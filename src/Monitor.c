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
  AppCtx*        user = (AppCtx*)ctx;   /* user-defined application context */
  PetscReal      dt;
  PetscReal      maxim, mean;
  PetscInt       timestep_number, size;

  // Compute field statistics:
  VecGetSize(u, &size);
  VecMax(u, NULL, &maxim);
  VecSum(u, &mean);
  mean /= size;
  TSGetTimeStep(ts,&dt);
  TSGetTimeStepNumber(ts,&timestep_number);

  // Abort if the solution is diverging
  if (maxim > 1e4 || maxim < -1e4) {
    PetscAbortErrorHandler(PETSC_COMM_WORLD, __LINE__ ,  __FUNCT__ , __FILE__ ,
                           PETSC_ERR_NOT_CONVERGED, 404,
                           "Solution is diverging over time.", NULL);
  }

  // Output to a *.vts file.
  if(user->param->write_each_step) {
    char filename[128];
    char* name = "output/solution";
    char num[5];
    char* ext  = ".vts";

    // Build the filename
    sprintf(num, "%1d", timestep_number);
    strncpy(filename, name, sizeof(filename));
    strncat(filename, num, (sizeof(filename) - strlen(filename)));
    strncat(filename, ext, (sizeof(filename) - strlen(filename)));

    PetscViewer viewer;
    PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
    VecView(u,viewer);
    PetscViewerDestroy(&viewer);
  }

  PetscPrintf(PETSC_COMM_WORLD,"step #: %D,\t step size = %.4f,\t time = %.2f,\t max = %.6f,\t mean = %.6f\n",step,(double)dt,(double)time, maxim, mean);

  return 0;
}
