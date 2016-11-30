/*
 * Monitor.c
 *
 *  Created on: Nov 26, 2016
 *      Author: clarkp
 */

#include "Monitor.h"
#include "Settings.h"

/* --------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Monitor"
/*
   Monitor - User-provided routine to monitor the solution computed at
   each timestep.  This example plots the solution and computes the
   error in two different norms.

   Input Parameters:
   ts     - the timestep context
   step   - the count of the current step (with 0 meaning the
             initial condition)
   time   - the current time
   u      - the solution at this timestep
   ctx    - the user-provided context for this monitoring routine.
            In this case we use the application context which contains
            information about the problem size, workspace and the exact
            solution.
*/
PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal time,Vec u, void *ctx) {
  AppCtx         *appctx = (AppCtx*)ctx;   /* user-defined application context */
  PetscReal      dt;
  PetscReal      maxim, mean;
  PetscInt       size;

  // Compute field statistics:
  VecGetSize(u, &size);
  VecMax(u, NULL, &maxim);
  VecSum(u, &mean);
  mean /= size;
  TSGetTimeStep(ts,&dt);

  // Abort if the solution is diverging
  if (maxim > 1e4 || maxim < -1e4) {
    PetscAbortErrorHandler(PETSC_COMM_WORLD, __LINE__ ,  __FUNCT__ , __FILE__ ,
                           PETSC_ERR_NOT_CONVERGED, 404,
                           "Solution is diverging over time.", NULL);
  }

  PetscPrintf(PETSC_COMM_WORLD,"step #: %D,\t step size = %g,\t time = %g,\t max = %g,\t mean = %g\n",step,(double)dt,(double)time, maxim, mean);

  return 0;
}
