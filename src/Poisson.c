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

// Set up the Poisson solver context only once, to avoid excessive memory
// allocation / deallocation
PetscErrorCode InitializePoissonContext(DM da, PoissonCtx* ctx, AppCtx* user) {
  KSPCreate(PETSC_COMM_WORLD, &ctx->ksp);

  // Set up the null space (used to make Neumann BCs consistent
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &ctx->nullspace);

  // Create the constant matrix for solving a Poisson equation
  DMCreateMatrix(da, &ctx->poisson_matrix);
  DMCreateMatrix(da, &ctx->poisson_matrix);
  SetUpPoissonMatrix(da, ctx->poisson_matrix, ctx->nullspace, user);
  KSPSetOperators(ctx->ksp, ctx->poisson_matrix, ctx->poisson_matrix);

  // Initialize the RHS matrix
  VecDuplicate(user->p, &ctx->RHS);

  KSPSetDM(ctx->ksp, da); // Sets the DM used by preconditioners
  KSPSetDMActive(ctx->ksp, PETSC_FALSE); // Tells KSP to keep our matrix
  // Use what's stored in the pressure as our initial guess:
  KSPSetInitialGuessNonzero(ctx->ksp,PETSC_TRUE);

  // Use the multigrid preconditioner and the BiCGSTAB solver
  PC pc;
  KSPGetPC(ctx->ksp, &pc);
  PCSetType(pc, PCGAMG);
  PCMGSetType(pc, PC_MG_KASKADE);
  KSPSetType(ctx->ksp, KSPBCGS);

  // Finish setting up ksp
  KSPSetFromOptions(ctx->ksp);
  KSPSetUp(ctx->ksp);

  return 0;
}

PetscErrorCode FreePoissonContext(PoissonCtx* ctx) {
  VecDestroy(&ctx->RHS);
  MatNullSpaceDestroy(&ctx->nullspace);
  MatDestroy(&ctx->poisson_matrix);
  KSPDestroy(&ctx->ksp);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SolvePoisson"
PetscErrorCode SolvePoisson(DM da_vel, DM da_p, PetscReal dt,
                              AppCtx *user,
                              PetscReal (*opt_func)(PetscReal, PetscReal)) {
  PetscLogEvent USER_EVENT;

  PetscLogEventRegister("Poisson Solver  ",0,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  // Create the RHS with the calculated divergence of u
  SetUpPoissonRHS(da_vel, da_p, dt, user->poisson_ctx->RHS,
                  user->poisson_ctx->nullspace, user, opt_func);

  // Solve the system
  KSPSolve(user->poisson_ctx->ksp, user->poisson_ctx->RHS, user->p);

  PetscLogEventEnd(USER_EVENT,0,0,0,0);

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

  hx = user->grid->dx;
  hy = user->grid->dy;
  mx = user->grid->mx;
  my = user->grid->my;

  // Get the start and size of the local array
  DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      // The row of the matrix corresponds to the i,j grid point
      row.i = i; row.j = j;
      if (i==0 || j==0 || i==mx-1 || j==my-1 ) {
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
        if (i!=mx-1) {
          v[num] = 1/(hx*hx);                col[num].i = i+1; col[num].j = j;
          numx++; num++;
        }
        if (j!=my-1) {
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
        v[4] = 1/(hy*hy);              col[4].i = i;   col[4].j = j+1;
        MatSetValuesStencil(poisson_matrix, 1, &row, 5, col, v, INSERT_VALUES);
      }
    }
  }

  // Fix the null space by setting one value to be = 0;
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
        dudx = (field[j][i+1].u - field[j][i].u)/hx;
        dvdy = (field[j+1][i].v - field[j][i].v)/hy;
        rhs[j][i] = (dudx + dvdy)/dt;
      } else {
        // NOTE: This section is intended for testing only.
        x = (i+.5)*hx;
        y = (j+.5)*hy;
        rhs[j][i] = opt_func(x, y);
      };
    }
  }

  // Put the array values back into the vectors
  DMDAVecRestoreArray    (da_p,   RHS,       &rhs);
  DMDAVecRestoreArrayRead(da_vel, local_vel, &field);

  // Normalize the output vector so the mean is zero
  PetscInt size;
  PetscReal sum;
  VecGetSize(RHS, &size);
  VecSum(RHS, &sum);
  VecShift(RHS, -1.0*sum/size);

  // Force the right hand side to be consistent for a singular matrix
  // Note that this is really a hack; normally the model would provide you with
  // a consistent RHS
  MatNullSpaceRemove(nullspace,RHS);

  DMRestoreLocalVector(da_vel, &local_vel);

  return 0;
}
