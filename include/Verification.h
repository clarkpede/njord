/*
 * Verification.h
 *
 *  Created on: Dec 3, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_VERIFICATION_H_
#define INCLUDE_VERIFICATION_H_

#include <petscdmda.h>
#include "Field.h"
#include "Settings.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Sets up a vector containing the exact pressure data.
 *
 * This function does NOT fill in the ghost points at the boundaries.
 *
 * --Inputs--
 * @param da - The distributed array matching the u and v velocities.
 * @param exact_u - A function pointer to a function specifying the analytical
 *                  solution to the u-velocity.
 * @param exact_v - A function pointer to a function specifying the analytical
 *                  solution to the v-velocity.
 * @param user - The user-specified context.
 *
 * --Outputs--
 * @param U - A global vector for the exact u and v velocities.
 * @return Returns 0 if the function finishes execution
 */
PetscErrorCode SetUpExactSolutionUV(DM da, Vec U, PetscReal time,
                                  PetscReal (*exact_u)(PetscReal, PetscReal, PetscReal),
                                  PetscReal (*exact_v)(PetscReal, PetscReal, PetscReal),
                                  AppCtx *user);

/** Sets up a vector containing the exact pressure data.
 *
 * This function does NOT fill in the ghost points at the boundaries.
 *
 * --Inputs--
 * @param da - The distributed array matching the pressure.
 * @param exact_p - A function pointer to a function specifying the analytical
 *                  solution to the pressure.
 * @param user - The user-specified context.
 *
 * --Outputs--
 * @param P - A global vector for the exact pressure data.
 * @return Returns 0 if the function finishes execution
 */
PetscErrorCode SetUpExactSolutionP(DM da, Vec P, PetscReal time,
                                  PetscReal (*exact_p)(PetscReal, PetscReal, PetscReal),
                                  AppCtx *user);

/** Computes the given error norm of the data.
 *
 * --Inputs--
 * @param exact - A vector representing the exact solution.
 * @param approx - A vector representing the approximate solution.
 * @param type - The norm type (NORM_1, NORM_2, or NORM_INFINITY)
 *
 * --Outputs--
 * @param err - The computed error norm.
 * @return Returns 0 if the function finishes execution.
 */
PetscErrorCode GetErrorNorm(Vec exact, Vec approx, NormType type,
                            PetscReal* err);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_VERIFICATION_H_ */
