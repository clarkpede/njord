/*
 * InterpolateInlet.h
 *
 *  Created on: Dec 6, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_INTERPOLATEINLET_H_
#define INCLUDE_INTERPOLATEINLET_H_

#include <petscsys.h>

#include "Settings.h"
#include "spline.h"

PetscErrorCode BuildInflowOutflowProfiles(PetscReal hy, PetscInt my,
                                          BoundaryProfiles* profiles);

PetscErrorCode FreeProfiles(BoundaryProfiles* profiles);

#endif /* INCLUDE_INTERPOLATEINLET_H_ */
