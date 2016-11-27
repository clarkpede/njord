/*
 * InitialSolution.h
 *
 *  Created on: Nov 23, 2016
 *      Author: clarkp
 */
#include <petscdmda.h>

#ifndef INCLUDE_INITIALSOLUTION_H_
#define INCLUDE_INITIALSOLUTION_H_

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode FormInitialFields(DM da, Vec x, PetscReal Lx, PetscReal Ly);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_INITIALSOLUTION_H_ */
