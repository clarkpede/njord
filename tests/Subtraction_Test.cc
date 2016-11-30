#include <boost/test/unit_test.hpp>
#include <petscsys.h>

#include <iostream>

BOOST_AUTO_TEST_SUITE(Subtraction)

BOOST_AUTO_TEST_CASE( a_minus_b ) {
  int a = 2;
  int b = 3;

  BOOST_CHECK(b-a == 1);
}

BOOST_AUTO_TEST_SUITE_END()
