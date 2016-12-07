/*
 * TimeMarch.c
 *
 *  Created on: Nov 23, 2016
 *      Author: clarkp
 *
 *-----------------------------------------------------------------------------
 * CONTROL VOLUMES
 *-----------------------------------------------------------------------------
 * U Control Volume
 *
   * v_i-1,j+1       v_i,j+1
   *  |^|--------------|^|
   *   |                |
   *   |                |
   * p_i-1,j   u_ij    p_ij
   *   |                |
   *   |                |
   *  |^|--------------|^|
   * v_i-1,j          v_ij
 *
 *-----------------------------------------------------------------------------
 *
 * V control volume
 *
 *   u_ij           u_i+1,j
 *    =>-----p_ij-----=>
 *    |                |
 *    |                |
 *    |      v_ij      |
 *    |                |
 *    |                |
 *    =>----p_i,j-1---=>
 *  u_i,j-1         u_i+1,j-1
 *
 *-----------------------------------------------------------------------------
 */

#include "TimeMarch.h"

PetscScalar get_FCu(Field **x, PetscReal hx, PetscReal hy,
                    PetscInt i, PetscInt j) {
  PetscReal      mN, mS, mE, mW;
  PetscReal      uN, uS, uE, uW;

    // Mass flow rates for u control volume
    mE =  0.5*hy*(x[j][i+1].u   + x[j][i].u);
    mW = -0.5*hy*(x[j][i-1].u   + x[j][i].u);
    mN =  0.5*hx*(x[j+1][i-1].v + x[j+1][i].v);
    mS = -0.5*hx*(x[j][i-1].v   + x[j][i].v);

    // Interpolated u velocities
    uE = 0.5*(x[j][i].u + x[j][i+1].u);
    uW = 0.5*(x[j][i].u + x[j][i-1].u);
    uN = 0.5*(x[j][i].u + x[j+1][i].u);
    uS = 0.5*(x[j][i].u + x[j-1][i].u);

  // U momentum convection, divided by area
  return  (mE*uE + mW*uW + mN*uN + mS*uS)/(hx*hy);
}


PetscScalar get_FCv(Field **x, PetscReal hx, PetscReal hy,
                    PetscInt i, PetscInt j) {
  PetscReal      mN, mS, mE, mW;
  PetscReal      vN, vS, vE, vW;

  // Mass flow rates for v control volume
  mE =  0.5*hy*(x[j-1][i+1].u   + x[j][i+1].u);
  mW = -0.5*hy*(x[j][i].u       + x[j-1][i].u);
  mN =  0.5*hx*(x[j][i].v       + x[j+1][i].v);
  mS = -0.5*hx*(x[j][i].v       + x[j-1][i].v);

  // Interpolated v velocities
  vE = 0.5*(x[j][i].v + x[j][i+1].v);
  vW = 0.5*(x[j][i].v + x[j][i-1].v);
  vN = 0.5*(x[j][i].v + x[j+1][i].v);
  vS = 0.5*(x[j][i].v + x[j-1][i].v);

  // V momentum convection, divided by area
  return (mE*vE + mW*vW + mN*vN + mS*vS)/(hx*hy);
}

PetscReal get_FDu(Field **x, PetscReal hx, PetscReal hy,
                  PetscInt i, PetscInt j, PetscReal nu) {
  PetscReal FDu;
  PetscReal dUdyN, dUdyS, dUdxE, dUdxW;

  // U diffusion from momentum equation:
  dUdxE = (x[j][i+1].u - x[j][i].u)/hx;
  dUdxW = (x[j][i].u   - x[j][i-1].u)/hx;
  dUdyN = (x[j+1][i].u - x[j][i].u)/hy;
  dUdyS = (x[j][i].u   - x[j-1][i].u)/hy;

  FDu = nu * (dUdxE*hy - dUdxW*hy + dUdyN*hx - dUdyS*hx);
  FDu /= (hx*hy);
  return FDu;
}

PetscReal get_FDv(Field **x, PetscReal hx, PetscReal hy,
                  PetscInt i, PetscInt j, PetscReal nu) {
  PetscReal FDv;
  PetscReal dVdyN, dVdyS, dVdxE, dVdxW;

  // V momentum equation:
  dVdxE = (x[j][i+1].v - x[j][i].v)/hx;
  dVdxW = (x[j][i].v   - x[j][i-1].v)/hx;
  dVdyN = (x[j+1][i].v - x[j][i].v)/hy;
  dVdyS = (x[j][i].v   - x[j-1][i].v)/hy;

  FDv = nu * (dVdxE*hy - dVdxW*hy + dVdyN*hx - dVdyS*hx);
  FDv /= (hx*hy);
  return FDv;
}

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
#define __FUNCT__ "FormTimeDerivativeFunction"
PetscErrorCode FormTimeDerivativeFunction(TS ts, PetscReal ftime, Vec X, Vec F, void *ptr) {
  AppCtx         *user = (AppCtx*)ptr;
  DM             da;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      hx,hy;
  PetscReal      FCu, FDu, FDv, FCv;
  Vec            localX;
  Field          **x, **f;
  enum BoundaryType bc;

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
      if (i==0 || i==Mx-1) {
        f[j][i].u = 0;
      } else {
        FCu =  get_FCu(x,hx,hy,i,j);
        FDu = get_FDu(x,hx,hy,i,j,user->param->nu);
        f[j][i].u = -FCu + FDu;
      }
      if (j==0 || j==My-1) {
        f[j][i].v = 0;
      } else {
        FCv = get_FCv(x,hx,hy,i,j);
        FDv = get_FDv(x,hx,hy,i,j,user->param->nu);
        f[j][i].v = -FCv + FDv;
      }
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
PetscErrorCode get_dt(DM da, Vec U, PetscReal* dt, AppCtx* user){
  PetscReal max_vel;
  PetscInt  mx, my;
  PetscReal hx, hy;
  PetscReal dt_conv, dt_diff;

  hx = user->grid->dx;

  // Stability limit for Euler forward with linear convection
  VecMax(U,NULL,&max_vel);
  dt_conv = user->param->CFL*hx/(max_vel);

  // Stability limit for Euler forward with linear diffusion
  dt_diff = user->param->CFL*(hx*hx)/user->param->nu;

  // Set the smaller of the two as an output
  *dt = ((dt_conv < dt_diff) ? dt_conv: dt_diff);

  return 0;
}

PetscErrorCode Prestep(TS ts) {
  PetscReal dt, time;
  DM        da_vel, da_p;
  AppCtx*   user;

  // Unpack the application state from TS
  TSGetApplicationContext(ts, &user);
  TSGetDM(ts, &da_vel);
  DMDAGetReducedDMDA(da_vel, 1, &da_p);
  DMSetApplicationContext(da_p,user);

  // Recompute the time step based on the CFL condition
  get_dt(da_vel, user->vel, &dt, user);
  TSSetTimeStep(ts, dt);

  // Set up profiling for the time-integration:
  PetscLogEventRegister("Time Integration",0,&user->current_event);
  PetscLogEventBegin(user->current_event,0,0,0,0);

  return 0;
}

PetscErrorCode Prestage(TS ts, PetscReal stagetime) {
  PetscReal dt, time;
  DM        da_vel, da_p;
  AppCtx*   user;

  // Unpack the application state from TS
  TSGetApplicationContext(ts, &user);
  TSGetDM(ts, &da_vel);
  DMDAGetReducedDMDA(da_vel, 1, &da_p);
  DMSetApplicationContext(da_p,user);

  TSGetTimeStep(ts, &dt);
  TSGetTime(ts,&time);

  UpdateBoundaryConditionsUV(da_vel, da_p, user->vel, user->p, time, dt, user);

  return 0;
}

PetscErrorCode PressureCorrection(TS ts) {
  PetscReal dt;
  DM da_vel, da_p;
  AppCtx* user;
  PetscInt timestep_number;

  //End profiling for the time-integration:
  TSGetApplicationContext(ts, &user);
  PetscLogEventEnd(user->current_event,0,0,0,0);

  // Unpack the application state from TS
  TSGetTimeStep(ts, &dt);
  TSGetTimeStepNumber(ts, &timestep_number);
  TSGetDM(ts, &da_vel);
  DMDAGetReducedDMDA(da_vel, 1, &da_p);
  DMSetApplicationContext(da_p,user);

  // Solves the Poisson equation and stores the result in user->p
  SolvePoisson(da_vel, da_p, dt, user, NULL);
  // Pressure ghost cells need to be updated to get grad(p) at left/bottom
  UpdateBoundaryConditionsP(da_p, user->p, user);
  CorrectVelocities(da_vel, da_p, dt, user);

  return 0;

};

/** Advances the solution forward in time.
 *
 * @param da - The distributed array used for the data.
 */
#undef __FUNCT__
#define __FUNCT__ "TimeMarch"
PetscErrorCode TimeMarch(TS* ts, DM da_vel, DM da_p, AppCtx *user) {
  SNES  snes;
  PetscReal dt;
  const PetscReal kStartTime = 0.0;
  const PetscReal kEndTime = 1.0;
  const PetscInt kMaxSteps = 30;

  TSType time_scheme;
  Mat Jac=NULL;
  Mat Jmf=NULL;

  //---------------------------------------------------------------------------
  // Set up timestepper
  //---------------------------------------------------------------------------
  TSCreate(PETSC_COMM_WORLD, ts);
  TSSetDM(*ts, da_vel);
  TSSetProblemType(*ts,TS_NONLINEAR);
  TSSetType(*ts, TSEULER);
  TSGetSNES(*ts,&snes);
  TSSetSolution(*ts, user->vel);

  TSSetRHSFunction(*ts,NULL,FormTimeDerivativeFunction,user);

  // Print information about the time-stepper
  TSMonitorSet(*ts,Monitor,user,NULL);

  // Set the initial time step
  get_dt(da_vel, user->vel, &dt, user);
  TSSetInitialTimeStep(*ts,kStartTime,dt);
  // The end time and max number of time-steps can also be set on the command
  // line using -ts_max_steps and -ts_final_time.
  TSSetDuration(*ts,kMaxSteps,kEndTime);
  TSSetExactFinalTime(*ts, TS_EXACTFINALTIME_MATCHSTEP);

  // Update with user-defined options
  TSSetApplicationContext(*ts, user);
  TSSetFromOptions(*ts);

  // Set up the Jacobian, if an implicit method is to be used
  TSGetType(*ts, &time_scheme);
  if (strcmp(time_scheme,TSCN) || strcmp(time_scheme,TSBEULER)) {
    DMCreateMatrix(da_vel,&Jac);
    MatCreateSNESMF(snes,&Jmf);
    SNESSetJacobian(snes,Jmf,Jac,SNESComputeJacobianDefaultColor,0);
  }

  // Print the information about the time-stepper
  if (user->param->verbose) TSView(*ts,PETSC_VIEWER_STDOUT_WORLD);

  // Use a pre-step so we can update BCs and set a variable time step
  // Use the pressure correction as a post-step
  TSSetPreStep(*ts, Prestep);
  TSSetPreStage(*ts,Prestage);
  TSSetPostStep(*ts, PressureCorrection);

  TSSolve(*ts, user->vel);

  return 0;
}

