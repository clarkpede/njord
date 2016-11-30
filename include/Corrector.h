/*
 * Corrector.h
 *
 *  Created on: Nov 29, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_CORRECTOR_H_
#define INCLUDE_CORRECTOR_H_

#include <petscdm.h>
#include <petscdmda.h>
#include "Settings.h"
#include "Field.h"

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode CorrectVelocities(DM da_vel, DM da_p, PetscReal dt,
                                 AppCtx *user);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_CORRECTOR_H_ */
