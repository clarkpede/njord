/*
 * Settings.c
 *
 *  Created on: Dec 1, 2016
 *      Author: clarkp
 */

#include "Settings.h"

PetscErrorCode SetGridDefaults(GridInfo *grid) {
  grid->mx = 240;
  grid->my = 60;
  grid->Lx = 36.0;
  grid->Ly = 9.0;
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
  param->write_each_step = PETSC_FALSE;
  strcpy(param->outfile,"output/average_fields.vts");
  return 0;
}
