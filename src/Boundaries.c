/*
 * Boundaries.c
 *
 *  Created on: Dec 1, 2016
 *      Author: clarkp
 */

#include "Boundaries.h"

PetscErrorCode GetBoundaryType(PetscInt i, PetscInt j,PetscInt mx,
                               PetscInt my, enum BoundaryType* type) {
  if (i==-1) {
    if (j==-1) {
      *type = BOTTOMLEFT;
    } else if (j==my) {
      *type = TOPLEFT;
    } else {
      *type = LEFT;
    }
  } else if (i==mx) {
    if (j==-1) {
      *type = BOTTOMRIGHT;
    } else if (j==my) {
      *type = TOPRIGHT;
    } else {
      *type = RIGHT;
    }
  } else if (j==-1) {
    *type = BOTTOM;
  } else if (j==my) {
    *type = TOP;
  } else {
    *type = NONE;
  }

  return 0;
}

PetscReal BC_U(PetscReal x, PetscReal y, PetscReal t, PetscReal nu) {
  return -cos(x)*sin(y)*exp(-2*nu*t);
}

PetscReal BC_V(PetscReal x, PetscReal y, PetscReal t, PetscReal nu) {
  return sin(x)*cos(y)*exp(-2*nu*t);
}

PetscErrorCode UpdateBoundaryConditionsUV(DM da_vel, DM da_p, Vec U, Vec P,
                                          PetscReal time, PetscReal dt,
                                          AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  Vec u_local, local_p;
  Field **field;
  PetscScalar **p;
  PetscReal hx, hy, x, y;
  PetscReal dpdx, dpdy;
  enum BoundaryType bc;

  DMDAGetInfo(da_vel, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
  // This is not the usual (Mx-1) because the staggered grid omits one of the
  // nodal points to maintain consistent u,v,p array sizes.
  hx = user->grid->Lx/(PetscReal)(mx);
  hy = user->grid->Ly/(PetscReal)(my);

  DMDAGetGhostCorners(da_vel,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries

  // Get a local array, including ghost points
  DMGetLocalVector(da_vel, &u_local);
  DMGlobalToLocalBegin(da_vel, U, INSERT_VALUES, u_local);
  DMGlobalToLocalEnd  (da_vel, U, INSERT_VALUES, u_local);
  DMDAVecGetArray(da_vel,u_local,&field);

  DMGetLocalVector(da_p, &local_p);
  DMGlobalToLocalBegin(da_p, P, INSERT_VALUES, local_p);
  DMGlobalToLocalEnd  (da_p, P, INSERT_VALUES, local_p);
  // The trailing "Read" on this function just mean *not* collective for MPI:
  DMDAVecGetArrayRead(da_p, local_p, &p);

  // x and y are locations of the center of the boundary, NOT the cell center
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      // Left BC and the corners have no effect on the solution
      GetBoundaryType(i,j,mx,my,&bc);
      if (bc == NONE) continue;

      // Calculate the gradient of the pressure-like term, for the correction
      // At the ghost points, indices lie off the main array.  We need +/-1
      if (bc == BOTTOM && i>0) {
        dpdx = (p[j+1][i] - p[j+1][i-1])/hx;
      } else if (bc == TOP && i>0) {
        dpdx = (p[j-1][i] - p[j-1][i-1])/hx;
      } else {
        dpdx = 0;
      }

      if (bc==BOTTOM) {
        x = i*hx; y = 0;
        field[j][i].u = 2*(BC_U(x,y,time,user->param->nu) - dt*dpdx)
                    - field[j+1][i].u;
      } else if (bc==TOP) {
        x = i*hx; y = my*hy;
        field[j][i].u = 2*(BC_U(x,y,time,user->param->nu) - dt*dpdx)
                    - field[j-1][i].u;
      } else if (bc==RIGHT) {
        // Used for pressure-poisson RHS calc
        x = mx*hx; y = (j+0.5)*hy;
        field[j][i].u = BC_U(x,y,time,user->param->nu) - dt*dpdx;
      }
    }
  }

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      GetBoundaryType(i,j,mx,my,&bc);
      // Corner and bottom BCs are not needed
      if (bc == NONE) continue;

      // Calculate the gradient of the pressure-like term, for the correction
      // At the ghost points, indices lie off the main array.  We need +/-1
      if (bc==LEFT && j>0) {
        dpdy = (p[j][i+1] - p[j-1][i+1])/hy;
      } else if (bc==RIGHT && j>0) {
        dpdy = (p[j][i-1] - p[j-1][i-1])/hy;
      } else {
        dpdy = 0;
      }

      if (bc==TOP) {
        // Used for pressure-poisson RHS calc
        y = my*hy; x = (i+0.5)*hx;
        field[j][i].v = (BC_V(x,y,time,user->param->nu) - dt*dpdy);
      } else if (bc==LEFT) {
        x = 0; y = j*hy;
        field[j][i].v = 2*(BC_V(x,y,time,user->param->nu) - dt*dpdy)
                - field[j][i+1].v;
      } else if (bc==RIGHT) {
        x = my*hy; y = j*hy;
        field[j][i].v = 2*(BC_V(x,y,time,user->param->nu) - dt*dpdy) - field[j][i-1].v;
      }
    }
  }

  DMDAVecRestoreArray(da_vel,u_local,&field);
  DMLocalToGlobalBegin(da_vel, u_local, INSERT_VALUES, U);
  DMLocalToGlobalEnd  (da_vel, u_local, INSERT_VALUES, U);
  DMRestoreLocalVector(da_vel, &u_local);

  DMDAVecRestoreArrayRead(da_p,   local_p,   &p);

  return 0;
};

PetscErrorCode UpdateBoundaryConditionsP(DM da_p, Vec P, AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  Vec p_local;
  PetscScalar **p;
  PetscReal hx, hy, x, y;
  enum BoundaryType bc;

  DMDAGetInfo(da_p, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
  // This is not the usual (Mx-1) because the staggered grid omits one of the
  // nodal points to maintain consistent u,v,p array sizes.
  hx = user->grid->Lx/(PetscReal)(mx);
  hy = user->grid->Ly/(PetscReal)(my);

  DMDAGetGhostCorners(da_p,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries

  // Get a local array, including ghost vectors
  DMGetLocalVector(da_p, &p_local);
  DMGlobalToLocalBegin(da_p, P, INSERT_VALUES, p_local);
  DMGlobalToLocalEnd  (da_p, P, INSERT_VALUES, p_local);
  DMDAVecGetArray(da_p,p_local,&p);

  //p is located at the center of each cell
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      // Corner, top, and right BCs aren't needed
      // Top and left aren't needed because of the missing cells
      // (No need to solve grad(p) for missing quantities)
      // All we're doing here is enforcing the Neumann BCs on pressure.
      GetBoundaryType(i,j,mx,my,&bc);
      if (bc==LEFT) {
        p[j][i] = p[j][i+1];
      } else if (bc==BOTTOM) {
        p[j][i] = p[j+1][i];
      }
    }
  }

  DMDAVecRestoreArray(da_p,p_local,&p);
  DMLocalToGlobalBegin(da_p, p_local, INSERT_VALUES, P);
  DMLocalToGlobalEnd  (da_p, p_local, INSERT_VALUES, P);
  DMRestoreLocalVector(da_p, &p_local);
  return 0;
}
