#include <boost/test/unit_test.hpp>
#include <petscsys.h>

BOOST_AUTO_TEST_SUITE(Addition)

BOOST_AUTO_TEST_CASE( a_plus_b )
{
  PetscMPIInt    rank, size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  std::cout << "Size: " << size << std::endl;

  int a = 2;
  int b = 3;

  BOOST_CHECK(a+b == 5);
}

BOOST_AUTO_TEST_SUITE_END()
