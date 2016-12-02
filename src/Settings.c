/*
 * Settings.c
 *
 *  Created on: Dec 1, 2016
 *      Author: clarkp
 */

#include "Settings.h"

PetscErrorCode SetGridDefaults(GridInfo *grid) {
  grid->mx = 64;
  grid->my = 64;
  grid->Lx = 2*M_PI;
  grid->Ly = 2*M_PI;
  // This is not the usual Lx/(mx-1) because we omit one set of points due
  // to the staggered grid arrangement.
  grid->dx = grid->Lx/grid->mx;
  grid->dy = grid->Ly/grid->my;
  grid->stencil = DMDA_STENCIL_STAR;
  grid->bc_x = DM_BOUNDARY_GHOSTED;
  grid->bc_y = DM_BOUNDARY_GHOSTED;
  grid->dof = 2;
  grid->stencil_width = 2;
  return 0;
}

PetscErrorCode SetParamDefaults(Parameters *param) {
  param->nu = 0.01;
  param->end_time = 1.0;
  param->CFL = 0.6;
  param->verbose = PETSC_FALSE;
  strcpy(param->outfile,"output/average_fields.vts");
  return 0;
}
