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
  grid->Lx = M_PI;
  grid->Ly = M_PI;
  grid->stencil = DMDA_STENCIL_STAR;
  grid->bc_x = DM_BOUNDARY_GHOSTED;
  grid->bc_y = DM_BOUNDARY_GHOSTED;
  grid->dof = 2;
  grid->stencil_width = 2;
  return 0;
}

PetscErrorCode SetParamDefaults(Parameters *param) {
  param->nu = 0.05;
  param->CFL = 0.6;
  param->verbose = PETSC_FALSE;
  param->write_output = PETSC_FALSE;
  strcpy(param->outfile,"output/average_fields.vts");
  return 0;
}
