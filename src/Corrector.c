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
  PetscLogEvent USER_EVENT;

  PetscLogEventRegister("Correct Vel.    ",0,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

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
      dpdx = (p[j][i]-p[j][i-1])/hx;
      dpdy = (p[j][i]-p[j-1][i])/hy;
      if (i != 0) { // Don't correct Dirichlet BC at inlet
        vel[j][i].u -= dt*dpdx;
      } else {
        vel[j][i].u = user->inlet_profile[j];
      }
      if (j != 0) { // Don't correct bottom wall
          vel[j][i].v -= dt*dpdy;
      } else {
        vel[j][i].v = 0;
      }
      // Top wall is not represented by our truncated grid
    }
  }

  // Put the array values back into the vectors
  DMDAVecRestoreArrayRead(da_p,   local_p,   &p);
  DMDAVecRestoreArray    (da_vel, user->vel, &vel);

  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return 0;
}
