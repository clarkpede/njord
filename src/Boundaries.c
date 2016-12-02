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

PetscErrorCode UpdateBoundaryConditionsUV(DM da_vel, Vec U, AppCtx *user) {
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

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      GetBoundaryType(i,j,mx,my,&bc);
      if (bc==BOTTOM) {
        field[j][i].u = -field[my-1][i].u;
        field[j][i].v = field[my-1][i].v;
      } else if (bc==TOP) {
        field[j][i].u = -field[0][i].u;
        field[j][i].v = field[0][i].v;
      } else if (bc==LEFT) {
        field[j][i].u = field[j][mx-1].u;
        field[j][i].v = -field[j][mx-1].v;
      } else if (bc==RIGHT) {
        field[j][i].u = field[j][0].u;
        field[j][i].v = -field[j][0].v;
      } else if (bc==TOPLEFT) {
        field[j][i].u = field[0][mx-1].u;
        field[j][i].v = field[0][mx-1].v;
      } else if (bc==BOTTOMLEFT) {
        field[j][i].u = field[my-1][mx-1].u;
        field[j][i].v = field[my-1][mx-1].v;
      } else if (bc==TOPRIGHT) {
        field[j][i].u = field[0][0].u;
        field[j][i].v = field[0][0].v;
      } else if (bc==BOTTOMRIGHT) {
        field[j][i].u = field[my-1][0].u;
        field[j][i].v = field[my-1][0].v;
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

  //vt is located at the center of each cell
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      GetBoundaryType(i,j,mx,my,&bc);
      if (bc==BOTTOM) {
        p[j][i] = p[my-1][i];
        p[j][i] = p[my-1][i];
      } else if (bc==TOP) {
        p[j][i] = p[0][i];
        p[j][i] = p[0][i];
      } else if (bc==LEFT) {
        p[j][i] = p[j][mx-1];
        p[j][i] = p[j][mx-1];
      } else if (bc==RIGHT) {
        p[j][i] = p[j][0];
        p[j][i] = p[j][0];
      }
    }
  }

  DMDAVecRestoreArray(da_p,p_local,&p);
  DMLocalToGlobalBegin(da_p, p_local, INSERT_VALUES, P);
  DMLocalToGlobalEnd  (da_p, p_local, INSERT_VALUES, P);
  return 0;
};
