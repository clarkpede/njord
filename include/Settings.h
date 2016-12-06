/*
 * Settings.h
 *
 *  Created on: Nov 25, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_SETTINGS_H_
#define INCLUDE_SETTINGS_H_

#include <petscdmda.h>

// All of the user-specified parameters that change execution
typedef struct {
  PetscReal nu;                        // Kinematic viscosity
  PetscReal CFL;                       // CFL number
  PetscBool verbose;                   // If true, print a full runtime summary
  PetscBool write_output;              // If true, output solution to *.vts
  char    outfile[PETSC_MAX_PATH_LEN]; // The name of the *.vts file to be used
} Parameters;

// The information on the grid
typedef struct {
  DMDAStencilType stencil;   // Stencil for distributed array (not FD)
  DMBoundaryType bc_x, bc_y; // Boundary type for distributed array (not FD)
  PetscInt mx, my;           // Grid resolution
  PetscReal Lx, Ly;          // Domain size
  PetscReal dx, dy;          // Grid spacing
  PetscInt stencil_width;
  PetscInt dof;              // For this problem 4: u, v, p, & vt
} GridInfo;

// The application context
typedef struct {
  Parameters  *param;    // Runtime parameters
  GridInfo    *grid;     // Grid information
  Vec         vel, p; // Vectors representing the solution at each time step
  PetscLogEvent current_event;
} AppCtx;

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode SetGridDefaults(GridInfo *grid);

PetscErrorCode SetParamDefaults(Parameters *param);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_SETTINGS_H_ */
