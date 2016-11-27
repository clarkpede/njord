#include "InitialSolution.h"

#define _USE_MATH_DEFINES // M_PI is the value of pi in <math.h>
#include <math.h>
#include "Field.h"

/** Creates the initial condition for the fields
 *
 * @param da
 * @param U
 * @param Lx
 * @param Ly
 * @return
 */
#undef __FUNCT__
#define __FUNCT__ "FormInitialFields"
PetscErrorCode FormInitialFields(DM da, Vec U, PetscReal Lx, PetscReal Ly) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  Field **u;
  PetscReal hx, hy, x, y;

  DMDAGetInfo(da, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
  // This is not the usual (Mx-1) because the staggered grid omits one of the
  // nodal points to maintain consistent u,v,p array sizes.
  hx = Lx/(PetscReal)(mx);
  hy = Ly/(PetscReal)(my);

  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries
  DMDAVecGetArray(da,U,&u);

  // p, vt are located at the center of each cell
  for (j=ys; j<ys+ym; j++) {
    y = j*hy + hy/2.0;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx + hx/2.0;
      u[j][i].p  = -0.25*(cos(2*x)+cos(2*y));
      u[j][i].vt = 0.0;
    }
  }

  // u is located at the east/west sides of each cell
  for (j=ys; j<ys+ym; j++) {
    y = j*hy + hy/2.0;
    for (i=xs; i<xs+xm; i++) {
      x = (i+1)*hx;
      u[j][i].u  = -cos(x)*sin(y);
    }
  }

  // v is located at the north/south sides of each cell
  for (j=ys; j<ys+ym; j++) {
    y = (j+1)*hy;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx + hx/2.0;
      u[j][i].v  = sin(x)*cos(y);
    }
  }

  DMDAVecRestoreArray(da,U,&u);
  return 0;
}
