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

typedef struct {
  PetscReal (*left)(PetscReal, PetscReal, PetscReal);
  PetscReal (*right)(PetscReal, PetscReal, PetscReal);
  PetscReal (*top)(PetscReal, PetscReal, PetscReal);
  PetscReal (*bottom)(PetscReal, PetscReal, PetscReal);
} BCs;

// Represents the type of boundary represented by a particular ghost cell.
enum BoundaryType {
  NONE,      //!< NONE
  TOP,       //!< TOP
  BOTTOM,    //!< BOTTOM
  LEFT,      //!< LEFT
  RIGHT,     //!< RIGHT
  TOPLEFT,   //!< TOPLEFT
  TOPRIGHT,  //!< TOPRIGHT
  BOTTOMLEFT,//!< BOTTOMLEFT
  BOTTOMRIGHT//!< BOTTOMRIGHT
};

/** Given the i,j index, this determines if a point belongs to a ghost cell and
 * if so, which boundary it lies on.
 *
 * --Inputs--
 * @param i - The x-index of the array
 * @param j - The y-index of the array
 * @param mx - The number of grid points in the x-direction.
 * @param my - The number of grid points in the y-direction.
 *
 * --Outputs--
 * @param type - An enum representing the type of boundary , or NONE if the
 *               point indicated by i,j is not a ghost cell.
 * @return Returns 0 if the function finishes execution.
 */
PetscErrorCode GetBoundaryType(PetscInt i, PetscInt j, PetscInt mx,
                               PetscInt my, enum BoundaryType* type);

/** Updates the ghost cells so that the boundary conditions on velocities
 *  are satisfied.
 *
 * --Inputs--
 * @param da_vel - The DM object used for the u and v velocities
 * @param U - The global vector representing u and v velocities
 * @param user - The user specified context
 *
 * --Outputs--
 * @return Returns 0 if the function finishes execution.
 */
PetscErrorCode UpdateBoundaryConditionsUV(DM da_vel, Vec U, PetscReal time,
                                          AppCtx *user);

/** Updates the ghost cells so that the homogeneous Neumann BCs on pressure
 * are satisfied.
 *
 * @param da_p - The DM object used for the pressure.
 * @param P - The global vector representing pressure.
 * @param user - The user specified context.
 * @return
 */
PetscErrorCode UpdateBoundaryConditionsP(DM da_p, Vec P, AppCtx *user);

PetscErrorCode SetUpTGVortexBCs(BCs* U_BCs, BCs* V_BCs);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_BOUNDARIES_H_ */
