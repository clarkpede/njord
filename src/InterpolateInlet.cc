/*
 * interpolate_inlet.cc
 *
 *  Created on: Dec 6, 2016
 *      Author: clarkp
 */

#include <petscsys.h>

#include "spline.h"

PetscErrorCode GetInflowU(PetscReal hy, PetscInt my, PetscReal *U) {
  // Hardcoded values from Driver and Seegmiller
  std::vector<double> y_data {1.00, 1.10, 1.15, 1.20, 1.30,
                              1.40, 1.50, 1.70, 2.00, 2.40,
                              2.80, 3.20, 3.60, 4.00, 5.00,
                              6.00, 7.50, 8.20, 9.00};
  std::vector<double> u_data {0.000, 0.700, 0.729, 0.757, 0.794,
                              0.825, 0.853, 0.895, 0.947, 0.998,
                              1.020, 1.019, 1.016, 1.014, 1.009,
                              1.007, 1.007, 0.940, 0.000};
  // Set up a cubic spline interpolant
  tk::spline interpolant;
  interpolant.set_points(y_data,u_data);

  // Build an interpolated array
  PetscReal y;
  for (PetscInt i=0; i<my; i++) {
    y = (i+0.5)*hy;
    if (y<1) {
      U[i] = 0.0;
    } else {
      U[i] = -4.0*(y-1.5)*(y-1.5) + 1;
      //U[i] = interpolant(y);
    }
  };

  return 0;
};


