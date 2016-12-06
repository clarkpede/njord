/*
 * RANSModels.h
 *
 *  Created on: Nov 20, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_RANSMODELS_H_
#define INCLUDE_RANSMODELS_H_

#include <petscdm.h>
#include <petscdmda.h>

class SAModel {
 public:
  SAModel();
  PetscReal get_residual(PetscScalar **u, PetscScalar **v,
                                  PetscScalar **nu_tilda, PetscInt i, PetscInt j,
                                  PetscReal d, PetscReal nu,
                                  PetscReal hx, PetscReal hy);
 private:
  PetscReal get_S(PetscScalar **u, PetscScalar **v,
                           PetscInt i, PetscInt j, PetscReal hx, PetscReal hy);
  // Model coefficients
  const double cw1;
  const double cw2;
  const double cw3;
  const double cb1;
  const double cb2;
  const double sigma;
  const double kappa;
  const double cv1;
  const double ct1;
  const double ct2;
  const double ct3;
  const double ct4;
};


#endif /* INCLUDE_RANSMODELS_H_ */
