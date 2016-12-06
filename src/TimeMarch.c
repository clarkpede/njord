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

PetscScalar get_FCu(Field **field, PetscReal hx, PetscReal hy, PetscInt i,
                    PetscInt j, PetscReal t, enum BoundaryType bc_type,
                    BCs U_BCs, BCs V_BCs) {
  PetscReal      mN, mS, mE, mW;
  PetscReal      uN, uS, uE, uW;
  PetscReal      x,y;

  // FIXME: Calculate convection for the right side
  // We don't need to calculate U convection for the left and right sides,
  // since U is prescribed there and U is located on those boundaries.
  if (bc_type == NONE || bc_type == TOP || bc_type == BOTTOM) {
    mW = -0.5*hy*(field[j][i-1].u   + field[j][i].u);
    mE =  0.5*hy*(field[j][i+1].u   + field[j][i].u);
    uW =  0.5*(field[j][i].u + field[j][i-1].u);
    uE =  0.5*(field[j][i].u + field[j][i+1].u);

    if (bc_type == TOP) {
      x = i*hx; y = (j+1)*hy;
      mN = hx*V_BCs.top(x,y,t);
      uN = U_BCs.top(x,y,t);
    } else {
      mN =  0.5*hx*(field[j+1][i-1].v + field[j+1][i].v);
      uN = 0.5*(field[j][i].u + field[j+1][i].u);
    }

    if (bc_type == BOTTOM) {
      x = i*hx; y = 0;
      mS = -hx*V_BCs.bottom(x,y,t);
      uS = U_BCs.bottom(x,y,t);
    } else {
      mS = -0.5*hx*(field[j][i-1].v   + field[j][i].v);
      uS = 0.5*(field[j][i].u + field[j-1][i].u);
    }

    return  (mE*uE + mW*uW + mN*uN + mS*uS)/(hx*hy);

  } else {
    // This branch won't normally be used, but just in case set FCu to 0.0.
    return 0.0;
  }
}


PetscScalar get_FCv(Field **field, PetscReal hx, PetscReal hy, PetscInt i,
                    PetscInt j, PetscReal t,
                    enum BoundaryType bc_type, BCs U_BCs, BCs V_BCs) {
  PetscReal      mN, mS, mE, mW;
  PetscReal      vN, vS, vE, vW;
  PetscReal      x,y;

  // FIXME: Calculate convection for the top side
  // We don't need to calculate V convection for the top and bottom,
  // since V is prescribed there and V is also located on those boundaries
  if (bc_type == NONE || bc_type == LEFT || bc_type == RIGHT) {
    mS = -0.5*hx*(field[j][i].v       + field[j-1][i].v);
    mN =  0.5*hx*(field[j][i].v       + field[j+1][i].v);
    vS =  0.5*(field[j][i].v + field[j-1][i].v);
    vN =  0.5*(field[j][i].v + field[j+1][i].v);

    if (bc_type==LEFT) {
      x = 0.0; y = j*hy;
      mW = -hy*U_BCs.left(x,y,t);
      vW = V_BCs.left(x,y,t);
    } else {
      mW = -0.5*hy*(field[j][i].u       + field[j-1][i].u);
      vW = 0.5*(field[j][i].v + field[j][i-1].v);
    }

    if (bc_type==RIGHT) {
      x = (i+1)*hx; y = j*hy;
      mE = hy*U_BCs.right(x,y,t);
      vE = V_BCs.right(x,y,t);
    } else {
      mE =  0.5*hy*(field[j-1][i+1].u   + field[j][i+1].u);
      vE = 0.5*(field[j][i].v + field[j][i+1].v);
    }

    return (mE*vE + mW*vW + mN*vN + mS*vS)/(hx*hy);

  } else {
    // This branch won't normally be used, but just in case set FCv to 0.0.
    return 0.0;
  }
}

PetscReal get_FDu(Field **field, PetscReal hx, PetscReal hy, PetscInt i,
                  PetscInt j, PetscReal nu,  PetscReal t,
                  enum BoundaryType bc_type, BCs U_BCs, BCs V_BCs) {
  PetscReal FDu;
  PetscReal dUdyN, dUdyS, dUdxE, dUdxW;
  PetscReal x,y;

  // We take one-sided derivatives for the top, bottom, and left.
  // NOTE: Most of these are 1st order accurate
  if (bc_type == TOP || bc_type == TOPRIGHT || bc_type == TOPLEFT) {
    x = i*hx; y = (j+1)*hy;
    dUdyN = (U_BCs.top(x,y,t) - field[j][i].u)/(0.5*hy);
  } else {
    dUdyN = (field[j+1][i].u - field[j][i].u)/hy;
  }

  if (bc_type == BOTTOM || bc_type == BOTTOMRIGHT || bc_type == BOTTOMLEFT) {
    x = i*hx; y = 0.0;
    dUdyS = (field[j][i].u - U_BCs.bottom(x,y,t))/(0.5*hy);
  } else {
    dUdyS = (field[j][i].u   - field[j-1][i].u)/hy;
  }

  if (bc_type == LEFT || bc_type == TOPLEFT || bc_type == BOTTOMLEFT) {
    // TODO: Use a higher order approximation here...
    dUdxW = (field[j][i+1].u - field[j][i].u)/(hx);
  } else {
    dUdxW = (field[j][i].u   - field[j][i-1].u)/hx;
  }

  // The right is special because of the missing nodes on the left.
  // 2nd order accurate
  if (bc_type == RIGHT || bc_type == TOPRIGHT || bc_type == BOTTOMRIGHT) {
    x = (i+1)*hx; y = (j+0.5)*hy;
    dUdxE = (U_BCs.right(x,y,t) - field[j][i].u)/(hx);
  } else {
    dUdxE = (field[j][i+1].u - field[j][i].u)/hx;
  }

  FDu = nu * (dUdxE*hy - dUdxW*hy + dUdyN*hx - dUdyS*hx);
  FDu /= (hx*hy);
  return FDu;
}

PetscReal get_FDv(Field **field, PetscReal hx, PetscReal hy, PetscInt i,
                  PetscInt j, PetscReal nu, PetscReal t,
                  enum BoundaryType bc_type, BCs U_BCs, BCs V_BCs) {
  PetscReal FDv;
  PetscReal dVdyN, dVdyS, dVdxE, dVdxW;
  PetscReal x,y;

  // We take one-sided derivatives for the top, bottom, and left.
  // NOTE: Most of these are 1st order accurate
  if (bc_type == BOTTOM || bc_type == BOTTOMRIGHT || bc_type == BOTTOMLEFT) {
    // TODO: Use a higher order approximation here...
    dVdyS = (field[j+1][i].v - field[j][i].v)/(hy);
  } else {
    dVdyS = (field[j][i].v   - field[j-1][i].v)/hy;
  }

  if (bc_type == LEFT || bc_type == TOPLEFT || bc_type == BOTTOMLEFT) {
    x = 0.0; y = j*hy;
    dVdxW = (field[j][i].v - V_BCs.left(x,y,t))/(0.5*hx);
  } else {
    dVdxW = (field[j][i].v   - field[j][i-1].v)/hx;
  }

  if (bc_type == RIGHT || bc_type == TOPRIGHT || bc_type == BOTTOMRIGHT) {
    x = (i+1)*hx; y = j*hy;
    dVdxE = (V_BCs.right(x,y,t) - field[j][i].v)/(0.5*hx);
  } else {
    dVdxE = (field[j][i+1].v - field[j][i].v)/hx;
  }

  // The top is special because of the missing v nodes.
  // 2nd order accurate
  if (bc_type == TOP || bc_type == TOPRIGHT || bc_type == TOPLEFT) {
    x = (i+0.5)*hx; y = (j+1)*hy;
    dVdyN = (V_BCs.top(x,y,t) - field[j][i].v)/(hy);
  } else {
    dVdyN = (field[j+1][i].v - field[j][i].v)/hy;
  }

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
  enum BoundaryType bc_type;
  BCs U_BCs, V_BCs;

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

  SetUpTGVortexBCs(&U_BCs, &V_BCs);

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
      GetBoundaryType(i,j,Mx,My,&bc_type);
      if (i==0 || i==Mx-1) {
        f[j][i].u = 0;
      } else {
        FCu =  get_FCu(x,hx,hy,i,j,ftime,bc_type,U_BCs,V_BCs);
        FDu = get_FDu(x,hx,hy,i,j,user->param->nu,ftime,bc_type,U_BCs,V_BCs);
        f[j][i].u = -FCu + FDu;
      }
      if (j==0 || j==My-1) {
        f[j][i].v = 0;
      } else {
        FCv = get_FCv(x,hx,hy,i,j,ftime,bc_type,U_BCs,V_BCs);
        FDv = get_FDv(x,hx,hy,i,j,user->param->nu,ftime,bc_type,U_BCs,V_BCs);
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

  *dt = dt_conv;

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

  TSGetTime(ts,&time);

  UpdateBoundaryConditionsUV(da_vel, user->vel, time, user);

  // Set up profiling for the time-integration:
  PetscLogEventRegister("Time Integration",0,&user->current_event);
  PetscLogEventBegin(user->current_event,0,0,0,0);

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

  // Output to a *.vts file.
  if(user->param->write_output) {
    char filename[128];
    char* name = "output/solution";
    char num[5];
    char* ext  = ".vts";

    // Build the filename
    sprintf(num, "%1d", timestep_number);
    strncpy(filename, name, sizeof(filename));
    strncat(filename, num, (sizeof(filename) - strlen(filename)));
    strncat(filename, ext, (sizeof(filename) - strlen(filename)));

    PetscViewer viewer;
    PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
    VecView(user->vel,viewer);
    PetscViewerDestroy(&viewer);
  }

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
  TSSetType(*ts, TSBEULER);
  TSGetSNES(*ts,&snes);
  TSSetSolution(*ts, user->vel);

  TSSetRHSFunction(*ts,NULL,FormTimeDerivativeFunction,user);

  // Print information about the time-stepper
  TSMonitorSet(*ts,Monitor,NULL,NULL);

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
  TSSetPostStep(*ts, PressureCorrection);

  TSSolve(*ts, user->vel);

  return 0;
}

