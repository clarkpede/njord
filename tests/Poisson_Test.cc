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

struct F{
  F() {
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
  ~F() {
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

PetscErrorCode SetUpExactSolution(DM da, Vec U, AppCtx *user) {
  PetscInt i, j, mx, my, xs, ys, xm, ym;
  PetscScalar **p;
  PetscReal hx, hy, x, y;

  hx = user->grid->dx;
  hy = user->grid->dy;

  DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL); // Get local grid boundaries
  DMDAVecGetArray(da,U,&p);

  for (j=ys; j<ys+ym; j++) {
    y = (j+.5)*hy;
    for (i=xs; i<xs+xm; i++) {
      x = (i+.5)*hx;
      p[j][i]  = -test_rhs_func(x,y);
    }
  }

  DMDAVecRestoreArray(da,U,&p);
  return 0;
};

BOOST_FIXTURE_TEST_SUITE(Poisson, F)

BOOST_AUTO_TEST_CASE( Poisson_Solution_of_Constant_Is_Zero ) {

  VecSet(user->vel,1.0);
  UpdateBoundaryConditionsUV(da_vel, user->vel, user);

  PetscReal arbitrary_dt = 1.0;
  SolvePoisson(da_vel, da_p, arbitrary_dt, user, NULL);

  PetscScalar tol = 1e-14;
  PetscReal sum;
  VecSum(user->p,&sum);
  BOOST_CHECK_SMALL(sum, tol);

}

BOOST_AUTO_TEST_CASE( Poisson_Solution_of_Trig_Is_Close ) {
  Vec exact;
  VecDuplicate(user->p, &exact);
  SetUpExactSolution(da_p, exact, user);

  //XXX: Print to file
  PetscViewer viewer;
  PetscViewerVTKOpen(PETSC_COMM_WORLD,"output/p_exact.vts",
                     FILE_MODE_WRITE,&viewer);
  VecView(exact,viewer);
  PetscViewerDestroy(&viewer);

  PetscReal arbitrary_dt = 1.0;
  SolvePoisson(da_vel, da_p, arbitrary_dt, user, &test_rhs_func);

  PetscScalar tol = 1e-3;
  PetscScalar norm;
  PetscInt size;

  //XXX: Print to file
  PetscViewerVTKOpen(PETSC_COMM_WORLD,"output/p.vts",
                     FILE_MODE_WRITE,&viewer);
  VecView(user->p,viewer);
  PetscViewerDestroy(&viewer);

  VecAXPY(user->p,-1.0,exact);
  VecNorm(user->p, NORM_2, &norm);
  VecGetSize(user->p, &size);
  BOOST_CHECK_SMALL(norm/size, tol);
}

BOOST_AUTO_TEST_SUITE_END()
