/*
 * Laplacian.h
 *
 *  Created on: Nov 28, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_LAPLACIAN_H_
#define INCLUDE_LAPLACIAN_H_

#include <petscdm.h>
#include <petscdmda.h>

#include "Settings.h"

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode SolveLaplacian(DM da_vel, DM da_p, PetscReal timestep,
                              AppCtx *user,
                              PetscReal (*opt_func)(PetscReal, PetscReal));

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_LAPLACIAN_H_ */
