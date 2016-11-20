#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Main_Runner"
#include <boost/test/unit_test.hpp>

#include <petscsys.h>

struct PETScConfig {
  PETScConfig() {
    MPI_Init(NULL, NULL);
    PetscInitialize(NULL, NULL, 0, 0);
  }
  ~PETScConfig() {
    PetscFinalize();
    MPI_Finalize();
  }
};

BOOST_GLOBAL_FIXTURE( PETScConfig )
