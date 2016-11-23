#include "RANSModels.h"

#include <algorithm>    // std::max and std::min
#include <cmath>

SAModel::SAModel() :  cb1(0.1355),
                      sigma(2.0/3),
                      cb2(0.622),
                      kappa(0.41),
                      cw2(0.3),
                      cw3(2.0),
                      cv1(7.1),
                      cl3(1.2),
                      cl4(0.5),
                      cw1(cb1/(kappa*kappa)+(1+cb2)/sigma) {
};

/**
 * Uses the SA model class to find the turbulent viscosity.
 *
 * @param nu_hat - The turbulent parameter.
 * @return The turbulent viscosity.
 */
double SAModel::get_nu_t(double nu_hat) {
  // FIXME: get nu_hat from advection equation;
  // FIXME: get nu from main scope
  return nu_hat*fv1;
};

/**
 * Sets the variables for the SA Model at a point
 *
 * The SA model involves both global constants and spatially
 * varying coefficients that are used in the transport equations.
 * In order to facillitate the use of these spatially varying constants,
 * all of them can be found (at a point) using this function.
 *
 * This function must be called for each new point being considered.
 *
 * @param nu_hat - The turbulent term
 * @param nu - The (molecular) kinematic viscosity
 * @param Omega - The magnitude of the vorticity.
 * @param dist - The distance to the nearest wall.
 */
void SAModel::set(double nu_hat, double nu, double Omega, double dist) {
  // FIXME: get nu_hat from advection equation;
  // FIXME: get nu from main scope
  d = dist;
  double chi = nu_hat/nu;
  fv1 = std::pow(chi,3)/(std::pow(chi,3)+std::pow(cv1,3));
  fv2 = 1.0 - chi/(1.0+chi*fv1);
  ft2 = cl3*std::exp(-cl4*chi*chi);
  // Spalart recommended in private communication limiting S to be no
  // smaller than 0.3*Omega
  S = std::max(Omega + nu_hat*fv2/(kappa*kappa*d*d), 0.3*Omega);
  r = std::min(nu_hat/(S*kappa*kappa*d*d),10.0);
  g = r + cw2*(std::pow(r,6) - r);
  fw = g*std::pow((1.0+std::pow(cw3,6))/
                  (std::pow(g,6)+std::pow(cw3,6)), (1.0/6));
};
