/*
 * InitialSolution.h
 *
 *  Created on: Nov 23, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_INITIALSOLUTION_H_
#define INCLUDE_INITIALSOLUTION_H_

#include <petscdmda.h>

#include "Settings.h"

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode SetInitialVelocities(DM da, Vec U, AppCtx *user);

PetscErrorCode SetInitialPressure(DM da, Vec U, AppCtx *user);

PetscErrorCode SetInitialVt(DM da, Vec U, AppCtx *user);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_INITIALSOLUTION_H_ */
