#include <math.h>
#include "petscksp.h"
#include "petscvec.h"
#include "petscmat.h"

#include "Laplacian.h"

/** This function is used to define the right-hand side of the Poisson equation
 *  to be solved
 */
double func(double x, double y) {
  return sin(x*M_PI)*sin(y*M_PI);
}

int main(int argc, char* argv[]) {
  KSP sles;
  Mat A;
  Vec b,x;
  int its, n;

  PetscInitialize(&argc, &argv, 0, 0);

  n = 10; // Mesh size is 10 by default
  PetscOptionsGetInt(PETSC_NULL, "-n", &n, 0);

  // Setup the problem
  A = FormLaplacian2d( n );
  b = FormVecFromFunction2d(n,func);
  VecDuplicate(b, &x);

  // Set up the solver
  KSPCreate(PETSC_COMM_WORLD, &sles);
  KSPSetOperators(sles, A, A);  // Should this include "DIFFERENT_NONZERO_PATTERN"?
  KSPSetFromOptions(sles);

  // Solve
  KSPSolve(sles, b, x);

  // Print results
  KSPGetIterationNumber(sles, &its);
  PetscPrintf(PETSC_COMM_WORLD, "Solution in %d iterations is \n",its);
  VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  // Clean up and exit
  MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
  KSPDestroy(&sles);
  PetscFinalize();
  return 0;

}
