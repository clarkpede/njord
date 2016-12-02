#include <math.h>
#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "Field.h"
#include "Settings.h"
#include "InitialSolution.h"
#include "TimeMarch.h"
#include "RANSModels.h"

static char help[] = "Solves 2D incompressible RANS equations using SA model.\n\n";

PetscErrorCode SetParams(Parameters*, GridInfo*);
PetscErrorCode ReportParams(Parameters*, GridInfo*);

// ----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char* argv[]) {
  Vec   x,r;      // Solution and residual vectors
  Mat   J;        // Jacobian matrix
  DM    da_vel, da_p;
  AppCtx *user;      // User defined work context
  Parameters param;  // Physical and other parameters
  GridInfo grid;     // Parameters defining the grid
  PetscErrorCode ierr;
  MPI_Comm comm;

  PetscInitialize(&argc, &argv,(char*)0, help);
  comm = PETSC_COMM_WORLD;

  // Set up the problem parameters
  SetParams(&param,&grid);
  ReportParams(&param,&grid);

  // Create user context and setup initial conditions
  PetscNew(&user);
  user->param = &param;
  user->grid =  &grid;

  // Create a distributed array for the velocities
  DMDACreate2d(comm, grid.bc_x, grid.bc_y, grid.stencil, grid.mx, grid.my,
               PETSC_DECIDE, PETSC_DECIDE, grid.dof, grid.stencil_width,
               0,0,&da_vel);
  DMSetApplicationContext(da_vel,user);
  DMDASetFieldName(da_vel,0,"x mean velocity");
  DMDASetFieldName(da_vel,1,"y mean velocity");
  DMCreateGlobalVector(da_vel,&(user->vel));

  // Create a distributed array for the pressure
  DMDAGetReducedDMDA(da_vel, 1, &da_p);
  DMSetApplicationContext(da_p,user);
  DMCreateGlobalVector(da_p,&(user->p));

  // Create the initial field
  SetInitialVelocities(da_vel, user->vel, user);
  SetInitialPressure(da_p, user->p, user);

  // Output the initial file to a *.vts file
  if(param.write_output) {
    PetscViewer viewer;
    PetscViewerVTKOpen(PETSC_COMM_WORLD,"output/initialfield.vts",
                       FILE_MODE_WRITE,&viewer);
    VecView(user->vel,viewer);
    PetscViewerDestroy(&viewer);
  }

  // Advance the solution through time
  TimeMarch(da_vel, da_p, user);

  // Output to a *.vts file.
  if(param.write_output) {
    PetscViewer viewer;
    PetscViewerVTKOpen(PETSC_COMM_WORLD,param.outfile,FILE_MODE_WRITE,&viewer);
    VecView(user->vel,viewer);
    PetscViewerDestroy(&viewer);
  }

  // Cleanup
  VecDestroy(&user->vel);
  VecDestroy(&user->p);
  PetscFree(user);
  DMDestroy(&da_vel);
  DMDestroy(&da_p);
  PetscPrintf(PETSC_COMM_WORLD,"Exited successfully.\n");
  PetscFinalize();
  return 0;
}
//-----------------------------------------------------------------------------


/** Sets up the runtime parameters using the options database
 *
 *  See:
 *  http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tutorials/ex30.c.html
 *  for an example of how to fill out this function.
 *
 * @param param - All runtime parameters will be stored in 'param'
 * @param grid - Grid specifications will be stored in 'grid'
 * @return Returns an error if execution fails
 */
PetscErrorCode SetParams(Parameters *param, GridInfo *grid) {
  PetscErrorCode ierr, ierr_out = 0;
  // Parameters
  PetscOptionsBegin(PETSC_COMM_WORLD,"","User-specified runtime options",
                    __FILE__);{
    param->nu = 0.0110;
    param->end_time = 1.0;
    PetscOptionsReal("-T","Duration of the simulation","none",param->end_time,
                     &(param->end_time),NULL);
    param->CFL = 0.6;
    PetscOptionsReal("-CFL","CFL number used for time advancement","none",
                     param->CFL,&(param->CFL),NULL);
    strcpy(param->outfile,"output/average_fields.vts");
    PetscOptionsString("-o","Output solution to a specified *.vts file",
                       "none",param->outfile,param->outfile, PETSC_MAX_PATH_LEN,
                       &param->write_output);
    param->verbose = PETSC_FALSE;
    PetscOptionsName("-verbose","Print all info during execution.","none",
                     &(param->verbose));
  }
  PetscOptionsEnd();

  // Grid Options
  PetscOptionsBegin(PETSC_COMM_WORLD,"","User-specified grid options",__FILE__); {
    grid->mx = 32;
    grid->my = 32;
    grid->Lx = 2*M_PI;
    grid->Ly = 2*M_PI;
    ierr = PetscOptionsInt("-mx","Grid resolution in x-direction","none",
                           grid->mx,&grid->mx,NULL);
    ierr = PetscOptionsInt("-my","Grid resolution in y-direction","none",
                           grid->my,&grid->my,NULL);
    // This is not the usual Lx/(mx-1) because we omit one set of points due
    // to the staggered grid arrangement.
    grid->dx = grid->Lx/grid->mx;
    grid->dy = grid->Ly/grid->my;
    grid->stencil = DMDA_STENCIL_STAR;
    grid->bc_x = DM_BOUNDARY_GHOSTED;
    grid->bc_y = DM_BOUNDARY_GHOSTED;
    grid->dof = 2;
    grid->stencil_width = 2;
  }
  PetscOptionsEnd();
  return ierr_out;

}

/** Prints out the relevant parameters, if full output is enabled.
 *
 * @param param - Runtime parameters used
 * @param grid - Grid information
 * @return Returns an error if execution fails
 */
PetscErrorCode ReportParams(Parameters *param, GridInfo *grid) {
  PetscErrorCode ierr, ierr_out=0;

  if (param->verbose) {
    PetscPrintf(PETSC_COMM_WORLD,
                "---------------------BEGIN NJORD PARAM REPORT-------------------\n");
    if (param->write_output) {
      PetscPrintf(PETSC_COMM_WORLD,"Writing final solution to: \n");
      PetscPrintf(PETSC_COMM_WORLD,"    %s", param->outfile);
    };
    PetscPrintf(PETSC_COMM_WORLD,"\nGrid: \n");
    PetscPrintf(PETSC_COMM_WORLD,"    [mx,my]: %D, %D \n", grid->mx, grid->my);
    PetscPrintf(PETSC_COMM_WORLD,
                "---------------------END NJORD PARAM REPORT-------------------\n");
  }
  return 0;
}
