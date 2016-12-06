#include "RANSModels.h"

#include <cmath>
#include <algorithm>

SAModel::SAModel() :  sigma(2.0/3),
                      cb1(0.1355),
                      cb2(0.622),
                      kappa(0.41),
                      cw1(cb1/(kappa*kappa)+(1+cb2)/sigma),
                      cw2(0.3),
                      cw3(2.0),
                      cv1(7.1),
                      ct1(1),
                      ct2(2),
                      ct3(1.2),
                      ct4(0.5) {
};


PetscReal SAModel::get_residual(PetscScalar **u, PetscScalar **v,
                                PetscScalar **nu_tilda, PetscInt i, PetscInt j,
                                PetscReal d, PetscReal nu,
                                PetscReal hx, PetscReal hy) {
 PetscReal S;
 PetscScalar chi, ft2, fv1, fv2, S_tilda, r, g, fw;
 PetscReal ddx, ddy, d2dx2, d2dy2;

 S = get_S(u,v,i,j,hx,hy);
 chi = nu_tilda[j][i]/nu;
 // I do not use the ft1 trip term
 ft2 = ct3*std::exp(-ct4*chi*chi);
 fv1 = std::pow(chi,3)/(std::pow(chi,3) + std::pow(cv1,3));
 fv2 = 1 - chi/(1+chi*fv1);
 // To avoid potential numerical problems, S_tilda should never be allowed
 // to reach zero or go negative.
 S_tilda = std::max(S + nu_tilda[j][i]/(kappa*kappa*d*d)*fv2,
                    0.3*S/std::sqrt(2));
 r = std::min(nu_tilda[j][i]/(S_tilda*kappa*kappa*d*d),10.0);
 g = r + cw2*(std::pow(r,6) - r);
 fw = g*std::pow((1+std::pow(cw3,6))/(std::pow(g,6)+std::pow(cw3,6)),1.0/6);

 // Calculate the 2nd order terms in nu_tilda
 // Viscosity is assumed to be constant, so \Del nu = 0
 ddx = (nu_tilda[j][i+1] - nu_tilda[j][i-1])/hx;
 ddy = (nu_tilda[j+1][i] - nu_tilda[j-1][i])/hy;
 d2dx2 = (nu_tilda[j][i+1] - 2*nu_tilda[j][i] + nu_tilda[j][i-1])/(hx*hx);
 d2dy2 = (nu_tilda[j][i+1] - 2*nu_tilda[j][i] + nu_tilda[j][i-1])/(hy*hy);

 // Build the residual
 PetscReal RHS = -(u[j][i] + u[j][i+1])*ddx - (v[j][i] + v[j+1][i])*ddy;
 RHS += cb1*(1-ft2)*S_tilda*nu_tilda[j][i];
 RHS += (nu + nu_tilda[j][i])*(d2dx2 + d2dy2)/sigma;
 RHS += (1+cb2)*(ddx*ddx + ddy*ddy)/sigma;
 RHS -= (cw1*fw - cb1/kappa*ft2)*std::pow(nu_tilda[j][i]/d,2);
 return 0.0;
}

PetscReal SAModel::get_S(PetscScalar **u, PetscScalar **v,
                         PetscInt i, PetscInt j, PetscReal hx, PetscReal hy) {
  PetscReal dudy, dvdx;
  PetscReal Omega_11, Omega_22, Omega_12, Omega_21;

  // Central differences
  dudy = ((u[j+1][i+1] + u[j+1][i]) - (u[j-1][i+1] + u[j-1][i]))/hy;
  dvdx = ((v[j+1][i+1] + v[j][i+1]) - (v[j+1][i-1] + v[j][i-1]))/hx;

  // Diagonal of the rotation matrix is zero.
  Omega_12 = 0.5*(dudy - dvdx);
  Omega_21 = -Omega_12;

  return std::sqrt(2.0*Omega_12*Omega_12 + 2.0*Omega_21*Omega_21);
}
