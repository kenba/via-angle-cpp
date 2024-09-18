//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024 Ken Barker. All Rights Reserved.
// (ken dot barker at via-technology dot co dot uk)
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// @file test_two_sum.cpp
/// @brief Contains tests for the two_sum function.
//////////////////////////////////////////////////////////////////////////////
#include "via/two_sum.hpp"
#include <boost/test/unit_test.hpp>

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_two_sum)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_two_sum_float) {
  // test double precision
  const auto [s1, t1] = via::two_sum(1.0F, -1e-08F);
  BOOST_CHECK_EQUAL(1.0F, s1);
  BOOST_CHECK_EQUAL(-1e-08F, t1);

  const auto [s2, t2] = via::two_sum(1.0F, 1.0F);
  BOOST_CHECK_EQUAL(2.0F, s2);
  BOOST_CHECK_EQUAL(0.0F, t2);

  const auto [s3, t3] = via::two_sum(0.2F, -1.0F);
  BOOST_CHECK_EQUAL(-0.80000000000000004F, s3);
  BOOST_CHECK_EQUAL(1.49011612e-08F, t3);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_two_sum_double) {
  // test double precision
  const auto [s1, t1] = via::two_sum(1.0, -1e-53);
  BOOST_CHECK_EQUAL(1.0, s1);
  BOOST_CHECK_EQUAL(-1e-53, t1);

  const auto [s2, t2] = via::two_sum(1.0, 1.0);
  BOOST_CHECK_EQUAL(2.0, s2);
  BOOST_CHECK_EQUAL(0.0, t2);

  const auto [s3, t3] = via::two_sum(0.2, -1.0);
  BOOST_CHECK_EQUAL(-0.80000000000000004, s3);
  BOOST_CHECK_EQUAL(5.5511151231257827e-17, t3);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_two_sum_long_double) {
  // test long double precision
  const auto [s1, t1] = via::two_sum(1.0L, -1e-53L);
  BOOST_CHECK_EQUAL(1.0L, s1);
  BOOST_CHECK_EQUAL(-1e-53L, t1);

  const auto [s2, t2] = via::two_sum(1.0L, 1.0L);
  BOOST_CHECK_EQUAL(2.0L, s2);
  BOOST_CHECK_EQUAL(0.0L, t2);

  const auto [s3, t3] = via::two_sum(0.2L, -1.0L);
  BOOST_CHECK_EQUAL(-0.800000000000000000011L, s3);
  BOOST_CHECK_EQUAL(1.35525271560688054251e-20L, t3);
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
