/*
 * Laplacian.cc
 *
 *  Created on: Nov 26, 2016
 *      Author: clarkp
 */

#include <petscksp.h>

#include "Laplacian.h"
#include "Field.h"

extern PetscErrorCode ComputeMatrix(KSP ksp, Mat J, Mat jac, void* ctx);
extern PetscErrorCode SetUpLaplacianRHS(DM da_vel, DM da_p, PetscReal dt,
                                        Vec RHS, AppCtx *user,
                                        PetscReal (*opt_func)(PetscReal, PetscReal));

#undef __FUNCT__
#define __FUNCT__ "SolveLaplacian"
PetscErrorCode SolveLaplacian(DM da_vel, DM da_p, PetscReal dt,
                              AppCtx *user,
                              PetscReal (*opt_func)(PetscReal, PetscReal)) {
  KSP ksp;
  Vec RHS;

  KSPCreate(PETSC_COMM_WORLD, &ksp);

  // Set up the linear algebra system
  VecDuplicate(user->p, &RHS);
  SetUpLaplacianRHS(da_vel, da_p, dt, RHS, user, opt_func);
  KSPSetComputeOperators(ksp, ComputeMatrix, user);

  // Set up the remaining KSP settings
  KSPSetDM(ksp, da_p);
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  // Solve the system
  KSPSolve(ksp,RHS,user->p);

  VecDestroy(&RHS);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SetUpLaplacianRHS"
PetscErrorCode SetUpLaplacianRHS(DM da_vel, DM da_p, PetscReal dt, Vec RHS,
                                 AppCtx *user,
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
        // NOTE: These x,y definitions are for testing only.
        // They do not match the rest of the staggered grid x,y definitions.
        x = i*hx;
        y = j*hy;
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
  MatNullSpace nullspace;
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
  MatNullSpaceRemove(nullspace,RHS);
  MatNullSpaceDestroy(&nullspace);

  return 0;
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
