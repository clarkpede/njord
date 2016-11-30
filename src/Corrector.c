/*
 * Corrector.c
 *
 *  Created on: Nov 29, 2016
 *      Author: clarkp
 */

#include "Corrector.h"

PetscErrorCode CorrectVelocities(DM da_vel, DM da_p, PetscReal dt,
                                 AppCtx *user) {

  Vec local_vel, local_p;
  Field **vel;
  PetscScalar **p;
  PetscReal hx, hy;
  PetscReal dpdx, dpdy;
  PetscInt xs,ys,xm,ym,i,j;

  hx = user->grid->dx;
  hy = user->grid->dy;

  // Get the local pressure vector and array, without ghost points
  DMGetLocalVector(da_p, &local_p);
  DMGlobalToLocalBegin(da_p, user->p, INSERT_VALUES, local_p);
  DMGlobalToLocalEnd  (da_p, user->p, INSERT_VALUES, local_p);
  // The trailing "Read" on this function just mean *not* collective for MPI:
  DMDAVecGetArrayRead(da_p, local_p, &p);

  // Set up the velocity with a local array, without ghost points
  DMDAVecGetArray(da_vel, user->vel, &vel);

  DMDAGetCorners(da_vel,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      // Add in the velocity correction.
      dpdx = (p[j][i+1]-p[j][i])/hx;
      dpdy = (p[j+1][i]-p[j][i])/hy;
      PetscPrintf(PETSC_COMM_WORLD,"dpdx:\t%g\t\tdpdy:\t%g\n", dpdx, dpdy);
      vel[j][i].u -= dt*dpdx;
      vel[j][i].v -= dt*dpdy;
    }
  }

  // Put the array values back into the vectors
  DMDAVecRestoreArrayRead(da_p,   local_p,   &p);
  DMDAVecRestoreArray    (da_vel, user->vel, &vel);

  return 0;
}
