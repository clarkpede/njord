/*
 * Laplacian.cc
 *
 *  Created on: Nov 26, 2016
 *      Author: clarkp
 */

#include <petscksp.h>

#include "Laplacian.h"

extern PetscErrorCode ComputeMatrix(KSP ksp, Mat J, Mat jac, void* ctx);
extern PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx);

#undef __FUNCT__
#define __FUNCT__ "SolveLaplacian"
PetscErrorCode SolveLaplacian(DM da_vel, DM da_p, PetscReal timestep,
                              AppCtx *user) {
  KSP ksp;
  Vec RHS;

  KSPCreate(PETSC_COMM_WORLD, &ksp);

  // Set up the linear algebra system
  KSPSetComputeRHS(ksp,ComputeRHS,user);
  KSPSetComputeOperators(ksp, ComputeMatrix, user);

  // Set up the reamining KSP settings
  KSPSetDM(ksp, da_p);
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  // Solve the system
  KSPSolve(ksp,NULL,NULL);
  KSPGetSolution(ksp,&(user->p));

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx)
{
  AppCtx    *user = (AppCtx*)ctx;
  PetscErrorCode ierr;
  PetscInt       i,j,mx,my,xm,ym,xs,ys;
  PetscScalar    Hx,Hy;
  PetscScalar    **rhs;
  DM             da;

  PetscFunctionBeginUser;
  KSPGetDM(ksp,&da);
  DMDAGetInfo(da, 0, &mx, &my, 0,0,0,0,0,0,0,0,0,0);
  Hx   = user->grid->Lx / (PetscReal)(mx);
  Hy   = user->grid->Ly / (PetscReal)(my);

  DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);
  DMDAVecGetArray(da, b, &rhs);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      rhs[j][i] = 0.0; // FIXME: How to handle u->p communication?
    }
  }
  DMDAVecRestoreArray(da, b, &rhs);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  /* force right hand side to be consistent for singular matrix */
  /* note this is really a hack, normally the model would provide you with a
   * consistent right-hand side */
  MatNullSpace nullspace;

  MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
  MatNullSpaceRemove(nullspace,b);
  MatNullSpaceDestroy(&nullspace);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(KSP ksp,Mat J,Mat jac,void *ctx)
{
  AppCtx         *user = (AppCtx*)ctx;
  PetscErrorCode ierr;
  PetscInt       i,j,mx,my,xm,ym,xs,ys;
  PetscScalar    v[5];
  PetscReal      Hx,Hy,HydHx,HxdHy,rho;
  MatStencil     row, col[5];
  DM             da;

  PetscFunctionBeginUser;
  ierr      = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  rho       = 1.0;
  ierr      = DMDAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  Hx        = user->grid->Lx / (PetscReal)(mx);
  Hy        = user->grid->Ly / (PetscReal)(my);
  HxdHy     = Hx/Hy;
  HydHx     = Hy/Hx;
  ierr      = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      row.i = i; row.j = j;
      if (i==0 || j==0 || i==mx|| j==my) {
        PetscInt numx = 0, numy = 0, num = 0;
        if (j!=0) {
          v[num] = -rho*HxdHy;              col[num].i = i;   col[num].j = j-1;
          numy++; num++;
        }
        if (i!=0) {
          v[num] = -rho*HydHx;              col[num].i = i-1; col[num].j = j;
          numx++; num++;
        }
        if (i!=mx) {
          v[num] = -rho*HydHx;              col[num].i = i+1; col[num].j = j;
          numx++; num++;
        }
        if (j!=my) {
          v[num] = -rho*HxdHy;              col[num].i = i;   col[num].j = j+1;
          numy++; num++;
        }
        v[num] = numx*rho*HydHx + numy*rho*HxdHy; col[num].i = i;   col[num].j = j;
        num++;
        ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);CHKERRQ(ierr);
      } else {
        v[0] = -rho*HxdHy;              col[0].i = i;   col[0].j = j-1;
        v[1] = -rho*HydHx;              col[1].i = i-1; col[1].j = j;
        v[2] = 2.0*rho*(HxdHy + HydHx); col[2].i = i;   col[2].j = j;
        v[3] = -rho*HydHx;              col[3].i = i+1; col[3].j = j;
        v[4] = -rho*HxdHy;              col[4].i = i;   col[4].j = j+1;
        ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // Fix the singular matrix problem created by the Neumann boundary conditions
  MatNullSpace nullspace;
  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
  ierr = MatSetNullSpace(J,nullspace);CHKERRQ(ierr);
  ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
