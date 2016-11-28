/*
 * TimeMarch.c
 *
 *  Created on: Nov 23, 2016
 *      Author: clarkp
 */
#include "math.h"

#include "TimeMarch.h"
#include "Field.h"
#include "Monitor.h"
#include "Laplacian.h"

/**
 *
 * @param ts - The TS context
 * @param ftime - The current time
 * @param X - The input vector, containing all state variables.
 * @param F - The function evaluation
 * @param ptr - Problem context, as set by SNESSetFunction()
 * @return Returns a zero if the function executes successfully.
 */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(TS ts, PetscReal ftime, Vec X, Vec F, void *ptr) {
  AppCtx         *user = (AppCtx*)ptr;
  DM             da;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      hx,hy;
  PetscReal      FCu, FDu, FDv, FCv;
  PetscReal      mN, mS, mE, mW;
  PetscReal      uN, uS, uE, uW;
  PetscReal      vN, vS, vE, vW;
  PetscReal      dVdyN, dVdyS, dVdxW, dVdxE;
  PetscReal      dUdxE, dUdxW, dUdyN, dUdyS;
  Vec            localX;
  Field          **x, **f;

  PetscFunctionBeginUser;
  TSGetDM(ts,&da);
  DMGetLocalVector(da,&localX);
  DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
              PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
              PETSC_IGNORE,PETSC_IGNORE);

  // This is not the usual (Mx-1) because the staggered grid omits one of the
  // nodal points to maintain consistent u,v,p array sizes.
  hx = user->grid->Lx/(PetscReal)(Mx);
  hy = user->grid->Ly/(PetscReal)(My);

  /*
       Scatter ghost points to local vector,using the 2-step process
          DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
       By placing code between these two statements, computations can be
       done while messages are in transition.
   */
  DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);
  DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);

  // Get pointers to vector data
  DMDAVecGetArrayRead(da,localX,&x);
  DMDAVecGetArray(da,F,&f);

  // Get local grid boundaries
  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);

  // Compute function over the locally owned part of the grid
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      //-----------------------------------------------------------------------
      // U advection:
      //-----------------------------------------------------------------------
      // Mass flow rates for u control volume
      mN =  0.5*hx*(x[j][i].v   + x[j][i+1].v);
      mS = -0.5*hx*(x[j-1][i].v + x[j-1][i+1].v);
      mE =  0.5*hy*(x[j][i].u   + x[j][i+1].u);
      mW = -0.5*hy*(x[j][i].u   + x[j][i-1].u);
      // Interpolated u velocities
      uN = 0.5*(x[j][i].u + x[j+1][i].u);
      uS = 0.5*(x[j][i].u + x[j-1][i].u);
      uE = 0.5*(x[j][i].u + x[j][i+1].u);
      uW = 0.5*(x[j][i].u + x[j][i-1].u);
      // U momentum convection, divided by area
      FCu =  (mN*uN + mS*uS + mE*uE + mW*uW)/(hx*hy);

      //-----------------------------------------------------------------------
      // V advection:
      //-----------------------------------------------------------------------
      // Mass flow rates for v control volume
      mN =  0.5*hx*(x[j][i].v   + x[j+1][i].v);
      mS = -0.5*hx*(x[j-1][i].v + x[j][i].v);
      mE =  0.5*hy*(x[j][i].u   + x[j+1][i].u);
      mW = -0.5*hy*(x[j][i-1].u + x[j+1][i-1].u);
      // Interpolated v velocities
      vN = 0.5*(x[j][i].v + x[j+1][i].v);
      vS = 0.5*(x[j][i].v + x[j-1][i].v);
      vE = 0.5*(x[j][i].v + x[j][i+1].v);
      vW = 0.5*(x[j][i].v + x[j][i-1].v);
      // V momentum convection, divided by area
      FCv =  (mN*vN + mS*vS + mE*vE + mW*vW)/(hx*hy);

      //------------------------------------------------------------------------
      // Diffusion, divided by the area
      //------------------------------------------------------------------------
      // U momentum equation:
      dUdyN = (x[j+1][i].u - x[j][i].u)/hy;
      dUdyS = (x[j][i].u   - x[j-1][i].u)/hy;
      dUdxE = (x[j][i+1].u - x[j][i].u)/hx;
      dUdxW = (x[j][i].u   - x[j][i-1].u)/hx;
      FDu = user->param->nu * (dUdyN*hx - dUdyS*hx + dUdxE*hy - dUdxW*hy);
      FDu /= (hx*hy);

      // V momentum equation:
      dVdyN = (x[j+1][i].v - x[j][i].v)/hy;
      dVdyS = (x[j][i].v   - x[j-1][i].v)/hy;
      dVdxE = (x[j][i+1].v - x[j][i].v)/hy;
      dVdxW = (x[j][i].v   - x[j][i-1].v)/hy;
      FDv = user->param->nu * (dVdyN*hx - dVdyS*hx + dVdxE*hy + dVdxE*hy);
      FDv /= (hx*hy);

      f[j][i].u = -FCu + FDu;
      f[j][i].v = -FCu + FDu;
    }
  }

  /*
       Restore vectors
   */
  DMDAVecRestoreArrayRead(da,localX,&x);
  DMDAVecRestoreArray(da,F,&f);
  DMRestoreLocalVector(da,&localX);
  PetscFunctionReturn(0);
}

/**
 *
 * @param da - The distributed array used for the data.
 * @param CFL - The CFL number used for each time step.
 * @param Lx - The length of the domain in the x-direction.
 * @param Ly - The length of the domain in the y-direction.
 * @param dt - A pointer to the time step to be set.
 * @return Returns a zero if the function executes successfully.
 */
#undef __FUNCT__
#define __FUNCT__ "get_dt"
PetscErrorCode get_dt(DM da, Vec U, PetscReal CFL,
                      PetscReal Lx, PetscReal Ly, PetscReal* dt){
  PetscReal max_vel;
  PetscInt  mx, my;
  PetscReal hx, hy;

  DMDAGetInfo(da, PETSC_IGNORE, &mx, &my, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
              PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);
  hx = Lx/(PetscReal)(mx);
  hy = Ly/(PetscReal)(my);

  VecMax(U,NULL,&max_vel);
  *dt = CFL*hx/(max_vel);

  return 0;
}

/** Advances the solution forward in time.
 *
 * @param da - The distributed array used for the data.
 */
#undef __FUNCT__
#define __FUNCT__ "TimeMarch"
PetscErrorCode TimeMarch(DM da_vel, DM da_p, DM da_vt, AppCtx *user) {
  TS    ts_conv, ts_diff;
  SNES  snes;
  PetscReal dt, time;
  PetscInt kMaxSteps = 1000;
  PetscInt n;


  //---------------------------------------------------------------------------
  // Set up timestepper
  //---------------------------------------------------------------------------
  TSCreate(PETSC_COMM_WORLD, &ts_conv);
  TSSetDM(ts_conv, da_vel);
  TSSetProblemType(ts_conv,TS_NONLINEAR);
  TSSetRHSFunction(ts_conv,NULL,FormFunction,user);
  TSSetType(ts_conv, TSEULER);
  TSGetSNES(ts_conv,&snes);
  TSSetSolution(ts_conv, user->vel);
  if (user->param->verbose) TSMonitorSet(ts_conv,Monitor,NULL,NULL);

  // Set the initial time step
  time = 0.0;
  get_dt(da_vel, user->vel, user->param->CFL,
         user->grid->Lx, user->grid->Ly, &dt);
  TSSetInitialTimeStep(ts_conv,time,dt);
  TSSetDuration(ts_conv,kMaxSteps,user->param->end_time);
  TSSetExactFinalTime(ts_conv, TS_EXACTFINALTIME_MATCHSTEP);

  // Update with user-defined options
  TSSetFromOptions(ts_conv);
  if (user->param->verbose) TSView(ts_conv,PETSC_VIEWER_STDOUT_SELF);

  n = 0;
  while (time < user->param->end_time) {
    if (time > 0.0) {
      // Recompute the time step based on the CFL condition
      get_dt(da_vel, user->vel, user->param->CFL,
             user->grid->Lx, user->grid->Ly, &dt);
      TSSetTimeStep(ts_conv, dt);
    }
    TSGetTimeStep(ts_conv,&dt); // Just in case dt is not set exactly

    Monitor(ts_conv,n,time,user->vel,NULL);
    TSStep(ts_conv);

    // Solve for the pressure-like term
    SolveLaplacian(da_vel, da_p, dt, user);

    // Correct the fields

    time += dt;
    n++;
  }

  return 0;
}

