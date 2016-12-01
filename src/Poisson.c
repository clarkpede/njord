/*
 * Poisson.c
 *
 *  Created on: Nov 26, 2016
 *      Author: clarkp
 */

#include <petscksp.h>

#include "Poisson.h"
#include "Field.h"

extern PetscErrorCode ComputeMatrix(KSP ksp, Mat J, Mat jac, void* ctx);
extern PetscErrorCode SetUpPoissonMatrix(DM da, Mat poisson_matrix,
                                         MatNullSpace nullspace, AppCtx *user);
extern PetscErrorCode SetUpPoissonRHS(DM da_vel, DM da_p, PetscReal dt,
                                        Vec RHS, MatNullSpace nullspace,
                                        AppCtx *user, PetscReal (*opt_func)
                                        (PetscReal, PetscReal));

#undef __FUNCT__
#define __FUNCT__ "SolvePoisson"
PetscErrorCode SolvePoisson(DM da_vel, DM da_p, PetscReal dt,
                              AppCtx *user,
                              PetscReal (*opt_func)(PetscReal, PetscReal)) {
  KSP ksp;
  Mat poisson_matrix;
  Vec RHS;
  MatNullSpace nullspace;

  KSPCreate(PETSC_COMM_WORLD, &ksp);

  // Set up the null space (used to make Neumann BCs consistent
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);

  // Create the constant matrix for solving a Poisson equation
  DMCreateMatrix(da_p, &poisson_matrix);
  SetUpPoissonMatrix(da_p, poisson_matrix, nullspace, user);
  KSPSetOperators(ksp, poisson_matrix, poisson_matrix);

  // Create the RHS with the calculated divergence of u
  VecDuplicate(user->p, &RHS);
  SetUpPoissonRHS(da_vel, da_p, dt, RHS, nullspace, user, opt_func);

  // Set up the remaining KSP settings
  KSPSetDM(ksp, da_p); // Sets the DM used by preconditioners
  KSPSetDMActive(ksp, PETSC_FALSE); // Tells KSP to keep our matrix
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  // Solve the system
  KSPSolve(ksp,RHS,user->p);

  // Normalize the output vector so the mean is zero
//  PetscInt size;
//  PetscReal sum;
//  VecGetSize(user->p, &size);
//  VecSum(user->p, &sum);
//  VecShift(user->p, -1.0*sum/size);

  VecDestroy(&RHS);
  MatNullSpaceDestroy(&nullspace);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SetUpPoissonMatrix"
PetscErrorCode SetUpPoissonMatrix(DM da, Mat poisson_matrix,
                                  MatNullSpace nullspace, AppCtx *user) {
  PetscInt       i,j,mx,my,xm,ym,xs,ys;
  PetscReal      hx, hy;
  MatStencil     row, col[5];
  PetscReal      v[5];

  hx = user->grid->dy;
  hy = user->grid->dx;

  // Get the start and size of the local array
  DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      // The row of the matrix corresponds to the i,j grid point
      row.i = i; row.j = j;
      if (i==0 || j==0 || i==mx|| j==my) {
        // We are on the boundary, and we need to use homogeneous Neumann BCs
        PetscInt numx = 0, numy = 0, num = 0;
        if (j!=0) {
          v[num] = 1/(hy*hy);                col[num].i = i;   col[num].j = j-1;
          numy++; num++;
        }
        if (i!=0) {
          v[num] = 1/(hx*hx);                col[num].i = i-1; col[num].j = j;
          numx++; num++;
        }
        if (i!=mx) {
          v[num] = 1/(hx*hx);                col[num].i = i+1; col[num].j = j;
          numx++; num++;
        }
        if (j!=my) {
          v[num] = 1/(hy*hy);                col[num].i = i;   col[num].j = j+1;
          numy++; num++;
        }
        v[num] = -numx/(hx*hx) - numy/(hy*hy); col[num].i = i;   col[num].j = j;
        num++;

        MatSetValuesStencil(poisson_matrix, 1, &row, num, col, v, INSERT_VALUES);

      } else {
        // We are inside the array, away from boundaries
        v[0] = 1/(hy*hy);              col[0].i = i;   col[0].j = j-1;
        v[1] = 1/(hx*hx);              col[1].i = i-1; col[1].j = j;
        v[2] = -2/(hx*hx) - 2/(hy*hy); col[2].i = i;   col[2].j = j;
        v[3] = 1/(hx*hx);              col[3].i = i+1; col[3].j = j;
        v[4] = 1/(hx*hx);              col[4].i = i;   col[4].j = j+1;
        MatSetValuesStencil(poisson_matrix, 1, &row, 5, col, v, INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(poisson_matrix,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(poisson_matrix,MAT_FINAL_ASSEMBLY);

  // Fix the singular matrix problem created by the Neumann boundary conditions
  MatSetNullSpace(poisson_matrix, nullspace);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SetUpPoissonRHS"
PetscErrorCode SetUpPoissonRHS(DM da_vel, DM da_p, PetscReal dt, Vec RHS,
                               MatNullSpace nullspace, AppCtx *user,
                               PetscReal (*opt_func)(PetscReal, PetscReal)){
  Vec local_vel, local_p;
  Field **field;
  PetscScalar **rhs;
  PetscReal hx, hy;
  PetscReal dudx, dvdy;
  PetscReal x,y;
  PetscInt xs,ys,xm,ym,i,j;

  hx = user->grid->dx;
  hy = user->grid->dy;

  // Get the local velocity vector and array, with ghost points
  DMGetLocalVector(da_vel, &local_vel);
  DMGlobalToLocalBegin(da_vel, user->vel, INSERT_VALUES, local_vel);
  DMGlobalToLocalEnd  (da_vel, user->vel, INSERT_VALUES, local_vel);
  // The trailing "Read" on this function just mean *not* collective for MPI:
  DMDAVecGetArrayRead (da_vel, local_vel, &field);

  // Set up the RHS with a local array, without ghost points
  DMDAVecGetArray(da_p, RHS, &rhs);

  DMDAGetCorners(da_p,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (opt_func == NULL) {
        dudx = (field[j][i].u - field[j][i-1].u)/hx;
        dvdy = (field[j][i].v - field[j-1][i].v)/hy;
        rhs[j][i] = (dudx + dvdy)/dt;
      } else {
        // NOTE: This section is intended for testing only.
        x = (i-.5)*hx;
        y = (j-.5)*hy;
        rhs[j][i] = opt_func(x, y);
      };
    }
  }

  // Put the array values back into the vectors
  DMDAVecRestoreArray    (da_p,   RHS,       &rhs);
  DMDAVecRestoreArrayRead(da_vel, local_vel, &field);

  // Force the right hand side to be consistent for a singular matrix
  // Note that this is really a hack; normally the model would provide you with
  // a consistent RHS
  MatNullSpaceRemove(nullspace,RHS);

  return 0;
}
