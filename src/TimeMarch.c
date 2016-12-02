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

      // Convection fluxes, divided by the area:
      FCu =  get_FCu(x,hx,hy,i,j);
      FCv =  get_FCv(x,hx,hy,i,j);

      // Diffusion fluxes, divided by the area:
      FDu = get_FDu(x,hx,hy,i,j,user->param->nu);
      FDv = get_FDv(x,hx,hy,i,j,user->param->nu);

      f[j][i].u = -FCu + FDu;
      f[j][i].v = -FCv + FDv;
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

  VecMax(U,NULL,&max_vel);
  dt_conv = user->param->CFL*hx/(max_vel);

  dt_diff = user->param->CFL*(hx*hx)/user->param->nu;

  // Return the smaller of the two
  *dt = ((dt_conv < dt_diff) ? dt_conv: dt_diff);

  return 0;
}

/** Advances the solution forward in time.
 *
 * @param da - The distributed array used for the data.
 */
#undef __FUNCT__
#define __FUNCT__ "TimeMarch"
PetscErrorCode TimeMarch(DM da_vel, DM da_p, AppCtx *user) {
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
  TSSetType(ts_conv, TSRK);
  TSGetSNES(ts_conv,&snes);
  TSSetSolution(ts_conv, user->vel);
  if (user->param->verbose) TSMonitorSet(ts_conv,Monitor,NULL,NULL);

  // Set the initial time step
  time = 0.0;
  get_dt(da_vel, user->vel, &dt, user);
  TSSetInitialTimeStep(ts_conv,time,dt);
  TSSetDuration(ts_conv,kMaxSteps,user->param->end_time);
  TSSetExactFinalTime(ts_conv, TS_EXACTFINALTIME_MATCHSTEP);

  // Update with user-defined options
  TSSetFromOptions(ts_conv);
  if (user->param->verbose) TSView(ts_conv,PETSC_VIEWER_STDOUT_WORLD);

  n = 0;
  while (time < user->param->end_time) {
    if (time > 0.0) {
      // Recompute the time step based on the CFL condition
      get_dt(da_vel, user->vel, &dt, user);
      TSSetTimeStep(ts_conv, dt);
    }
    TSGetTimeStep(ts_conv,&dt); // In case dt is not set exactly, read it

    UpdateBoundaryConditionsUV(da_vel, user->vel, user);
    UpdateBoundaryConditionsP(da_p, user->p, user);

    Monitor(ts_conv,n,time,user->vel,NULL); // Manually print out monitor data
    TSStep(ts_conv);

    // Solves the Poisson equation and stores the result in user->p
    SolvePoisson(da_vel, da_p, dt, user, NULL);
    CorrectVelocities(da_vel, da_p, dt, user);

    // Output to a *.vts file.
    if(user->param->write_output) {
      char filename[128];
      char* name = "output/solution";
      char num[5];
      char* ext  = ".vts";

      // Build the filename
      sprintf(num, "%1d", n);
      strncpy(filename, name, sizeof(filename));
      strncat(filename, num, (sizeof(filename) - strlen(filename)));
      strncat(filename, ext, (sizeof(filename) - strlen(filename)));

      PetscViewer viewer;
      PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
      VecView(user->vel,viewer);
      PetscViewerDestroy(&viewer);
    }

    time += dt;
    n++;
  }

  return 0;
}

