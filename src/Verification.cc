/*
 * Verification.cc
 *
 *  Created on: Dec 3, 2016
 *      Author: clarkp
 */

#include "Verification.h"

PetscErrorCode SetUpExactSolutionUV(DM da, Vec U,
                                  PetscReal (*exact_u)(PetscReal, PetscReal),
                                  PetscReal (*exact_v)(PetscReal, PetscReal),
                                  AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  Field **field;
  PetscReal hx, hy, x, y;

  hx = user->grid->dx;
  hy = user->grid->dy;

  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries
  DMDAVecGetArray(da,U,&field);

  for (j=ys; j<ys+ym; j++) {
    y = (j+.5)*hy;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx;
      field[j][i].u = (*exact_u)(x,y);
    }
  }

  for (j=ys; j<ys+ym; j++) {
    y = j*hy;
    for (i=xs; i<xs+xm; i++) {
      x = (i+.5)*hx;
      field[j][i].v = (*exact_v)(x,y);
    }
  }

  DMDAVecRestoreArray(da,U,&field);
  return 0;
};


PetscErrorCode SetUpExactSolutionP(DM da, Vec P,
                                  PetscReal (*exact_p)(PetscReal, PetscReal),
                                  AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  PetscReal **p;
  PetscReal hx, hy, x, y;

  hx = user->grid->dx;
  hy = user->grid->dy;

  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries
  DMDAVecGetArray(da,P,&p);

  for (j=ys; j<ys+ym; j++) {
    y = (j+.5)*hy;
    for (i=xs; i<xs+xm; i++) {
      x = (i+.5)*hx;
      p[j][i]  = (*exact_p)(x,y);
    }
  }

  DMDAVecRestoreArray(da,P,&p);
  return 0;
};


PetscErrorCode GetErrorNorm(Vec exact, Vec approx, NormType type,
                            PetscReal* err) {
  PetscScalar norm;
  PetscInt size;

  VecAXPY(approx,-1.0,exact);
  VecNorm(approx, type, &norm);
  VecGetSize(approx, &size);
  *err = norm/size;

  return 0;
};


