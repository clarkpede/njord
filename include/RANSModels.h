/*
 * RANSModels.h
 *
 *  Created on: Nov 20, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_RANSMODELS_H_
#define INCLUDE_RANSMODELS_H_

class SAModel {
 public:
  SAModel();
  void set(double nu_hat, double nu, double Omega, double dist);
  double get_nu_t(double nu_hat);
  const double cw1;
  const double cw2;
  const double cw3;
  const double cb1;
  const double cb2;
  const double sigma;
  const double kappa;
  const double cv1;
  const double cl3;
  const double cl4;
  double d;
  double fv1;
  double fv2;
  double fw;
  double g;
  double r;
  double ft2;
  double S;
};


#endif /* INCLUDE_RANSMODELS_H_ */
