/*
 * Monitor.h
 *
 *  Created on: Nov 26, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_MONITOR_H_
#define INCLUDE_MONITOR_H_

#include <petscdmda.h>
#include <petscts.h>

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal time,Vec u,void *ctx);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_MONITOR_H_ */
