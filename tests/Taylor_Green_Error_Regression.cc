#include <boost/test/unit_test.hpp>
#include <petscsys.h>

#include <petscdm.h>
#include <petscdmda.h>

#include <cmath>

#include "Settings.h"
#include "InitialSolution.h"
#include "TimeMarch.h"
#include "Verification.h"

BOOST_AUTO_TEST_SUITE(Taylor_Green_Error_Regression)

PetscReal exact_u(PetscReal x, PetscReal y, PetscReal t) {
  return -std::cos(x)*std::sin(y)*std::exp(-2*t);
}

PetscReal exact_v(PetscReal x, PetscReal y, PetscReal t) {
  return std::sin(x)*std::cos(y)*std::exp(-2*t);
}

void RunTest(PetscInt mx, PetscInt my, PetscReal* L2_err, PetscReal* Linf_err) {
  DM    da_vel, da_p;
  TS    ts;
  AppCtx *user;      // User defined work context
  Parameters param;  // Physical and other parameters
  GridInfo grid;     // Parameters defining the grid

  // Set up the problem parameters
  SetParamDefaults(&param);
  SetGridDefaults(&grid);
  grid.mx = mx; grid.my = my;
  param.nu = 1.0;
  param.CFL = 1.0;
  grid.dx = grid.Lx/grid.mx;
  grid.dy = grid.Ly/grid.my;

  // Create user context and setup initial conditions
  PetscNew(&user);
  user->param = &param;
  user->grid =  &grid;

  // Create a distributed array for the velocities
  DMDACreate2d(PETSC_COMM_WORLD, grid.bc_x, grid.bc_y, grid.stencil,
               grid.mx, grid.my, PETSC_DECIDE, PETSC_DECIDE,
               grid.dof, grid.stencil_width, 0,0,&da_vel);
  DMSetApplicationContext(da_vel,user);
  DMDASetFieldName(da_vel,0,"x mean velocity");
  DMDASetFieldName(da_vel,1,"y mean velocity");
  DMCreateGlobalVector(da_vel,&(user->vel));

  // Create a distributed array for the pressure
  DMDAGetReducedDMDA(da_vel, 1, &da_p);
  DMSetApplicationContext(da_p,user);
  DMCreateGlobalVector(da_p,&(user->p));

  // Create the initial field
  SetInitialVelocities(da_vel, user->vel, user);
  SetInitialPressure(da_p, user->p, user);

  // Advance the solution through time
  // NOTE: The number of steps is fixed at 30 here, to match Kim & Moin's paper
  PetscOptionsSetValue(NULL,"-ts_final_time","10000");
  PetscOptionsSetValue(NULL,"-ts_max_steps","30");
  TimeMarch(&ts, da_vel, da_p, user);

  // Compute the Linf error
  Vec exact;
  PetscReal time;
  TSGetTime(ts,&time);
  VecDuplicate(user->vel, &exact);
  SetUpExactSolutionUV(da_vel, exact, time, &exact_u, &exact_v, user);
  GetErrorNorm(exact, user->vel, NORM_2, L2_err);
  GetErrorNorm(exact, user->vel, NORM_INFINITY, Linf_err);
  PetscPrintf(PETSC_COMM_WORLD,"Linf Error:\t%g\n",*Linf_err);

  // Cleanup
  VecDestroy(&user->vel);
  VecDestroy(&user->p);
  PetscFree(user);
  DMDestroy(&da_vel);
  DMDestroy(&da_p);
}

BOOST_AUTO_TEST_CASE( Taylor_Green_Spatial_Convergence) {
  PetscReal L2_err[3], Linf_err[3];

  RunTest(16,16,&L2_err[0],&Linf_err[0]);
  RunTest(32,32,&L2_err[1],&Linf_err[1]);
  RunTest(64,64,&L2_err[2],&Linf_err[2]);

  PetscReal L2_rate_1 = std::log(L2_err[0]/L2_err[1])/std::log(2);
  PetscReal L2_rate_2 = std::log(L2_err[1]/L2_err[2])/std::log(2);
  PetscPrintf(PETSC_COMM_WORLD,"L2 Convergence rates:\t%g and %g\n",
              L2_rate_1, L2_rate_2);

  PetscReal Linf_rate_1 = std::log(Linf_err[0]/Linf_err[1])/std::log(2);
  PetscReal Linf_rate_2 = std::log(Linf_err[1]/Linf_err[2])/std::log(2);
  PetscPrintf(PETSC_COMM_WORLD,"Linf Convergence rates:\t%g and %g\n",
              Linf_rate_1, Linf_rate_2);

  double tol = .05;
  BOOST_CHECK_SMALL(L2_rate_1-2.0,tol);
  BOOST_CHECK_SMALL(L2_rate_2-2.0,tol);
  BOOST_CHECK_SMALL(Linf_rate_1-2.0,tol);
  BOOST_CHECK_SMALL(Linf_rate_2-2.0,tol);
}

BOOST_AUTO_TEST_SUITE_END()
