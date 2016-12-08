#include "InitialSolution.h"

#define _USE_MATH_DEFINES // M_PI is the value of pi in <math.h>
#include <math.h>
#include "Field.h"

/** Creates the initial condition for the fields
 *
 * @param da
 * @param U
 * @param user
 * @return
 */
#undef __FUNCT__
#define __FUNCT__ "SetInitialVelocities"
PetscErrorCode SetInitialVelocities(DM da, Vec U, AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  Field **u;
  PetscReal hx, hy, x, y;
  PetscReal blending;
  PetscReal outlet_profile;

  DMDAGetInfo(da, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
  // This is not the usual (Mx-1) because the staggered grid omits one of the
  // nodal points to maintain consistent u,v,p array sizes.
  hx = user->grid->Lx/(PetscReal)(mx);
  hy = user->grid->Ly/(PetscReal)(my);

  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries
  DMDAVecGetArray(da,U,&u);

  // u is located at the east/west sides of each cell
  for (j=ys; j<ys+ym; j++) {
    y = j*hy + hy/2.0;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx;
      blending = pow((1.0-x/user->grid->Lx),5.0);
      u[j][i].u  =         blending*user->profiles->inlet_u[j]
                   + (1.0-blending)*user->profiles->outlet_u[j];
    }
  }

  // v is located at the north/south sides of each cell
  for (j=ys; j<ys+ym; j++) {
    y = j*hy;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx + hx/2.0;
      blending = pow((1.0-x/user->grid->Lx),5.0);
      u[j][i].v  =         blending*user->profiles->inlet_v[j]
                   + (1.0-blending)*user->profiles->outlet_v[j];
    }
  }

  DMDAVecRestoreArray(da,U,&u);
  return 0;
}

/** Creates the initial condition for the fields
 *
 * @param da
 * @param U
 * @param user
 * @return
 */
#undef __FUNCT__
#define __FUNCT__ "SetInitialPressure"
PetscErrorCode SetInitialPressure(DM da, Vec U, AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  PetscScalar **p;
  PetscReal hx, hy, x, y;

  DMDAGetInfo(da, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
  // This is not the usual (Mx-1) because the staggered grid omits one of the
  // nodal points to maintain consistent u,v,p array sizes.
  hx = user->grid->Lx/(PetscReal)(mx);
  hy = user->grid->Ly/(PetscReal)(my);

  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries
  DMDAVecGetArray(da,U,&p);

  // p is located at the center of each cell
  for (j=ys; j<ys+ym; j++) {
    y = j*hy + hy/2.0;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx + hx/2.0;
      p[j][i]  = 0.0;
    }
  }

  DMDAVecRestoreArray(da,U,&p);
  return 0;
}

/** Creates the initial condition for the fields
 *
 * @param da
 * @param U
 * @param user
 * @return
 */
#undef __FUNCT__
#define __FUNCT__ "SetInitialPressure"
PetscErrorCode SetInitialVt(DM da, Vec U, AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  PetscScalar **vt;
  PetscReal hx, hy, x, y;

  DMDAGetInfo(da, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
  // This is not the usual (Mx-1) because the staggered grid omits one of the
  // nodal points to maintain consistent u,v,p array sizes.
  hx = user->grid->Lx/(PetscReal)(mx);
  hy = user->grid->Ly/(PetscReal)(my);

  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries
  DMDAVecGetArray(da,U,&vt);

  //vt is located at the center of each cell
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      vt[j][i] = 0.0;
    }
  }

  DMDAVecRestoreArray(da,U,&vt);
  return 0;
}
