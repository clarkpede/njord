#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <petscsys.h>

#define __USE_MATH_DEFINES
#include <cmath>

#include <petscdm.h>
#include <petscdmda.h>

#include "Poisson.h"
#include "Boundaries.h"
#include "Settings.h"
#include "Verification.h"

struct poisson_fixture{
  poisson_fixture() {
    //--------------------------
    // Setup
    //--------------------------
    // Build the Settings struct
    PetscNew(&user);
    user->grid = &grid;
    user->param = &params;

    SetGridDefaults(&grid);
    SetParamDefaults(&params);
    user->grid->dx = user->grid->Lx/user->grid->mx;
    user->grid->dy = user->grid->Ly/user->grid->my;

    // Build the DA objects
    // Create a distributed array for the velocities
    DMDACreate2d(PETSC_COMM_WORLD, grid.bc_x, grid.bc_y, grid.stencil, grid.mx, grid.my,
                 PETSC_DECIDE, PETSC_DECIDE, grid.dof, grid.stencil_width,
                 0,0,&da_vel);
    DMSetApplicationContext(da_vel,user);
    DMDASetFieldName(da_vel,0,"x mean velocity");
    DMDASetFieldName(da_vel,1,"y mean velocity");
    DMCreateGlobalVector(da_vel,&(user->vel));

    // Create a distributed array for the pressure
    DMDAGetReducedDMDA(da_vel, 1, &da_p);
    DMSetApplicationContext(da_p,user);
    DMCreateGlobalVector(da_p,&(user->p));
  }
  //--------------------------
  // Teardown
  //--------------------------
  ~poisson_fixture() {
    VecDestroy(&user->vel);
    VecDestroy(&user->p);
    PetscFree(user);
    DMDestroy(&da_vel);
    DMDestroy(&da_p);
  }

  Parameters params;
  GridInfo grid;
  AppCtx* user;
  DM da_vel;
  DM da_p;

};

PetscReal test_rhs_func(PetscReal x, PetscReal y) {
  // Test function has to satisfy Neumann BCs
  return std::cos(x) + std::cos(y);
};

PetscReal exact_solution(PetscReal x, PetscReal y, PetscReal t) {
  return -std::cos(x) - std::cos(y);
};

BOOST_FIXTURE_TEST_SUITE(Poisson, poisson_fixture)

BOOST_AUTO_TEST_CASE( Poisson_Solution_of_Constant_Is_Zero ) {

  VecSet(user->vel,0.0);
  // Ghost cells are automatically zero, so we don't need to update BCs

  PetscReal arbitrary_dt = 1.0;
  SolvePoisson(da_vel, da_p, arbitrary_dt, user, NULL);

  PetscScalar tol = 1e-14;
  PetscReal sum;
  VecSum(user->p,&sum);
  BOOST_CHECK_SMALL(sum, tol);

}

BOOST_AUTO_TEST_CASE( Poisson_Solution_of_Trig_Is_Close ) {
  PetscReal arbitrary_time = 0.0;

  Vec exact;
  VecDuplicate(user->p, &exact);
  SetUpExactSolutionP(da_p, exact, arbitrary_time, &exact_solution, user);

  PetscReal arbitrary_dt = 1.0;
  SolvePoisson(da_vel, da_p, arbitrary_dt, user, &test_rhs_func);

  PetscReal err;
  GetErrorNorm(exact, user->p, NORM_2, &err);

  PetscScalar tol = 1e-5;
  BOOST_CHECK_SMALL(err, tol);
}

BOOST_AUTO_TEST_SUITE_END()
