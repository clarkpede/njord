/*
 * TimeMarch.h
 *
 *  Created on: Nov 23, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_TIMEMARCH_H_
#define INCLUDE_TIMEMARCH_H_

#include <petscdmda.h>
#include <petscts.h>

#include "math.h"

#include "Poisson.h"
#include "Corrector.h"
#include "Field.h"
#include "Monitor.h"
#include "Settings.h"

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode TimeMarch(DM da_vel, DM da_p, AppCtx *user);

#ifdef __cplusplus
}
#endif



#endif /* INCLUDE_TIMEMARCH_H_ */
