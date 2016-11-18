/*
 * Laplacian.c
 *
 *  Created on: Nov 18, 2016
 *      Author: clarkp
 */

#include "Laplacian.h"

Mat FormLaplacian2d(int n) {
  Mat A;
  int r, rowStart, rowEnd, i, j;
  double h, oneByh2;

  h = 1.0/(n+1);
  oneByh2 = 1.0/(h*h);

  MatCreate(PETSC_COMM_WORLD, &A);

  MatSetSizes(A,PETSC_DECIDE, PETSC_DECIDE, n*n, n*n);
  // Set basic options (such as data structure) from the command line:
  MatSetFromOptions(A);
  MatSetUp(A);
  MatGetOwnershipRange(A, &rowStart, &rowEnd);

  for (r=rowStart; r<rowEnd; r++) {
    i = r%n; j = r/n;
    // Braces are needed here because of MatSetValue macro expansion...
    if (j-1 > 0) {
      MatSetValue(A,r,r-n,oneByh2, INSERT_VALUES);
    }
    if (i-1 > 0) {
      MatSetValue(A,r,r-1,oneByh2, INSERT_VALUES);
    }
    MatSetValue(A,r,r,-4*oneByh2,INSERT_VALUES);
    if (i+1<n-1) {
      MatSetValue(A,r,r+1,oneByh2,INSERT_VALUES);
    }
    if (j+1<n-1) {
      MatSetValue(A,r,r+n,oneByh2,INSERT_VALUES);
    }
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  return A;
}

Vec FormVecFromFunction2d(int n, double (*f)(double, double)){
  Vec V;
  int r, rowStart, rowEnd, i, j;
  double h;

  h = 1.0/(n+1);

  VecCreate(PETSC_COMM_WORLD,&V);
  VecSetSizes(V, PETSC_DECIDE, n*n);
  VecSetFromOptions(V);

  VecGetOwnershipRange(V, &rowStart, &rowEnd);
  for (r=rowStart; r<rowEnd; r++) {
    i = (r%n)+1;
    j = (r/n)+1;
    VecSetValue(V,r,(*f)(i*h, j*h),INSERT_VALUES);
  }

  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
  return V;

}


