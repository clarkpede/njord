/*
 * Boundaries.h
 *
 *  Created on: Dec 1, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_BOUNDARIES_H_
#define INCLUDE_BOUNDARIES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <petscdmda.h>

#include "math.h"

#include "Settings.h"
#include "Field.h"

enum BoundaryType {
  NONE,
  TOP,
  BOTTOM,
  LEFT,
  RIGHT,
  TOPLEFT,
  TOPRIGHT,
  BOTTOMLEFT,
  BOTTOMRIGHT
};

PetscErrorCode GetBoundaryType(PetscInt i, PetscInt j, PetscInt mx,
                               PetscInt my, enum BoundaryType* type);

PetscErrorCode UpdateBoundaryConditionsUV(DM da_vel, Vec U, AppCtx *user);

PetscErrorCode UpdateBoundaryConditionsP(DM da_p, Vec P, AppCtx *user);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_BOUNDARIES_H_ */
