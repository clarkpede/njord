/*
 * Boundaries.c
 *
 *  Created on: Dec 1, 2016
 *      Author: clarkp
 */

#include "Boundaries.h"

PetscErrorCode GetBoundaryType(PetscInt i, PetscInt j,PetscInt mx,
                               PetscInt my, enum BoundaryType* type) {
  if (i==0) {
    if (j==0) {
      *type = BOTTOMLEFT;
    } else if (j==my-1) {
      *type = TOPLEFT;
    } else {
      *type = LEFT;
    }
  } else if (i==mx-1) {
    if (j==0) {
      *type = BOTTOMRIGHT;
    } else if (j==my-1) {
      *type = TOPRIGHT;
    } else {
      *type = RIGHT;
    }
  } else if (j==0) {
    *type = BOTTOM;
  } else if (j==my-1) {
    *type = TOP;
  } else {
    *type = NONE;
  }

  return 0;
}

PetscReal BC_U(PetscReal x, PetscReal y, PetscReal t) {
  PetscReal nu = 1.0;
  return -cos(x)*sin(y)*exp(-2*nu*t);
}

PetscReal BC_V(PetscReal x, PetscReal y, PetscReal t) {
  PetscReal nu = 1.0;
  return sin(x)*cos(y)*exp(-2*nu*t);
}

PetscErrorCode SetUpTGVortexBCs(BCs* U_BCs, BCs* V_BCs) {
  PetscReal nu = 1;

  U_BCs->left = BC_U;
  U_BCs->right = BC_U;
  U_BCs->top = BC_U;
  U_BCs->bottom = BC_U;

  V_BCs->left = BC_V;
  V_BCs->right = BC_V;
  V_BCs->top = BC_V;
  V_BCs->bottom = BC_V;

  return 0;
};



PetscErrorCode UpdateBoundaryConditionsUV(DM da_vel, Vec U, PetscReal time,
                                          AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  Vec u_local;
  Field **field;
  PetscReal hx, hy, x, y;
  enum BoundaryType bc;

  DMDAGetInfo(da_vel, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
  // This is not the usual (Mx-1) because the staggered grid omits one of the
  // nodal points to maintain consistent u,v,p array sizes.
  hx = user->grid->Lx/(PetscReal)(mx);
  hy = user->grid->Ly/(PetscReal)(my);

  DMDAGetGhostCorners(da_vel,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries

  // Get a local array, including ghost vectors
  DMGetLocalVector(da_vel, &u_local);
  DMGlobalToLocalBegin(da_vel, U, INSERT_VALUES, u_local);
  DMGlobalToLocalEnd  (da_vel, U, INSERT_VALUES, u_local);
  DMDAVecGetArray(da_vel,u_local,&field);

//  // x and y are locations of the center of the boundary, NOT the cell center
//  for (j=ys; j<ys+ym; j++) {
//    for (i=xs; i<xs+xm; i++) {
//      // Left BC and the corners have no effect on the solution
//      GetBoundaryType(i,j,mx,my,&bc);
//      if (bc == NONE) continue;
//      if (bc==BOTTOM) {
//        x = i*hx; y = 0;
//        field[j][i].u = 2*BC_U(x,y,time) - field[j+1][i].u;
//      } else if (bc==TOP) {
//        x = i*hx; y = my*hy;
//        field[j][i].u = 2*BC_U(x,y,time) - field[j-1][i].u;
//      } else if (bc==RIGHT) {
//        // Used for pressure-poisson RHS calc
//        x = mx*hx; y = (j+0.5)*hy;
//        field[j][i].u = BC_U(x,y,time);
//      }
//    }
//  }

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      GetBoundaryType(i,j,mx,my,&bc);
      // Corner and bottom BCs are not needed
      if (bc == NONE) continue;
      if (bc==TOP) {
        // Used for pressure-poisson RHS calc
        y = my*hy; x = (i+0.5)*hx;
        field[j][i].v = BC_V(x,y,time);
      } else if (bc==LEFT) {
        x = 0; y = j*hy;
        field[j][i].v = 2*BC_V(x,y,time) - field[j][i+1].v;
      } else if (bc==RIGHT) {
        x = my*hy; y = j*hy;
        field[j][i].v = BC_V(x,y,time);
      }
    }
  }

  DMDAVecRestoreArray(da_vel,u_local,&field);
  DMLocalToGlobalBegin(da_vel, u_local, INSERT_VALUES, U);
  DMLocalToGlobalEnd  (da_vel, u_local, INSERT_VALUES, U);
  DMRestoreLocalVector(da_vel, &u_local);
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
