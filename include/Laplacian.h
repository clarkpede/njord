/*
 * Laplacian.h
 *
 *  Created on: Nov 18, 2016
 *      Author: clarkp
 */

#ifndef INCLUDE_LAPLACIAN_H_
#define INCLUDE_LAPLACIAN_H_

#include "petscvec.h"
#include "petscmat.h"

#ifdef __cplusplus
extern "C" {
#endif


Mat FormLaplacian2d(int n);

Vec FormVecFromFunction2d(int n, double (*f)(double, double));


#ifdef __cplusplus
}
#endif

#endif /* INCLUDE_LAPLACIAN_H_ */
