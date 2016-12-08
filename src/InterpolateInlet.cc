/*
 * interpolate_inlet.cc
 *
 *  Created on: Dec 6, 2016
 *      Author: clarkp
 */

#include "InterpolateInlet.h"

PetscErrorCode InterpolateProfiles(PetscReal hy, PetscInt my,
                                   BoundaryProfiles *profiles) {
  // Hardcoded values from Driver and Seegmiller
  std::vector<double> y_init {
    1.000, 1.100, 1.150, 1.200, 1.300, 1.400,
    1.500, 1.700, 2.000, 2.400, 2.800, 3.200,
    3.600, 4.000, 5.000, 6.000, 7.500, 8.200,
    9.000};

  std::vector<double> u_init {
    0.000, 0.700, 0.729, 0.757, 0.794, 0.825,
    0.853, 0.895, 0.947, 0.998, 1.020, 1.019,
    1.016, 1.014, 1.009, 1.007, 1.007, 0.940,
    0.000};

  std::vector<double> v_init {
    0.000, -0.017, -0.017, -0.020, -0.020, -0.020,
    -0.022, -0.023, -0.022, -0.020, -0.017, -0.016,
    -0.016, -0.017, -0.022, -0.020, -0.022, -0.015,
    0.000};

  std::vector<double> y_exit {
    0.000, 0.100, 0.150, 0.200, 0.300, 0.400,
    0.500, 0.600, 0.700, 0.800, 0.900, 1.000,
    1.100, 1.150, 1.200, 1.300, 1.400, 1.500,
    1.700, 2.000, 2.400, 2.800, 3.200, 3.600,
    4.000, 5.000, 6.000, 7.500, 9.000};

  std::vector<double> u_exit {
    0.000, 0.581, 0.604, 0.621, 0.645, 0.662,
    0.679, 0.691, 0.705, 0.718, 0.734, 0.747,
    0.761, 0.769, 0.778, 0.792, 0.808, 0.821,
    0.848, 0.886, 0.922, 0.936, 0.941, 0.942,
    0.943, 0.944, 0.944, 0.912, 0.000};

  std::vector<double> v_exit {
    0.000, -0.005, -0.004, -0.004, -0.005, -0.004,
    -0.005, -0.006, -0.006, -0.008, -0.009, -0.010,
    -0.009, -0.011, -0.013, -0.013, -0.013, -0.013,
    -0.012, -0.011, -0.008, -0.006, -0.003, -0.004,
    -0.005, -0.006, -0.007, -0.005, 0.000};


  // Set up a cubic spline interpolant
  tk::spline inlet_u_interpolant, inlet_v_interpolant;
  tk::spline outlet_u_interpolant, outlet_v_interpolant;

  // Set up the interpolation
  inlet_u_interpolant.set_points(y_init, u_init);
  inlet_v_interpolant.set_points(y_init, v_init);
  inlet_u_interpolant.set_points(y_exit, u_exit);
  inlet_u_interpolant.set_points(y_exit, v_exit);

  // Build interpolated arrays for the inlet
  PetscReal y;
  PetscInt  j;
  profiles->total_flux_in = 0;
  for (j=0; j<my; j++) {
    y = (j+0.5)*hy;
    if (y<1) {
      profiles->inlet_u[j] = 0.0;
      profiles->inlet_v[j] = 0.0;
    } else {
      //profiles->inlet_u[j] = inlet_u_interpolant(y);
      profiles->inlet_u[j] = -4.0*(y-1.5)*(y-1.5) + 1;
      profiles->inlet_v[j] = 0.0;
      profiles->total_flux_in += profiles->inlet_u[j];
    }
  };

  // Build interpolated arrays for the outlet
  for (j=0; j<my; j++) {
    y = (j+0.5)*hy;
    //profiles->inlet_u[j] = inlet_u_interpolant(y);
    profiles->outlet_u[j] = -0.5*(y-1)*(y-1) + 0.5;
    profiles->outlet_v[j] = 0.0;
  };

  return 0;
};

PetscErrorCode BuildInflowOutflowProfiles(PetscReal hy, PetscInt my,
                                          BoundaryProfiles* profiles) {
  PetscMalloc1(my, &profiles->inlet_u);
  PetscMalloc1(my, &profiles->inlet_v);
  PetscMalloc1(my, &profiles->outlet_u);
  PetscMalloc1(my, &profiles->outlet_v);

  InterpolateProfiles(hy, my, profiles);

  return 0;
};

PetscErrorCode FreeProfiles(BoundaryProfiles* profiles) {
  PetscFree(profiles->inlet_u);
  PetscFree(profiles->inlet_v);
  PetscFree(profiles->outlet_u);
  PetscFree(profiles->outlet_v);

  return 0;
}




