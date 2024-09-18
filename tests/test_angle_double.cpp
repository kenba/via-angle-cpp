//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024 Ken Barker. All Rights Reserved.
// (ken dot barker at via-technology dot co dot uk)
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// @file test_angle_doub;e.cpp
/// @brief Contains unit tests for the via angle classes and functions.
//////////////////////////////////////////////////////////////////////////////
#include "via/angle.hpp"
#include <boost/test/unit_test.hpp>

using namespace via;

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_angle)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_degrees) {
  const auto one{Degrees(1.0)};
  const auto two{Degrees(2.0)};
  const auto m_two{-two};

  BOOST_CHECK_EQUAL(-2.0, m_two.v());
  BOOST_CHECK_EQUAL(two.v(), m_two.abs().v());

  const auto m_one{one + m_two};
  BOOST_CHECK_EQUAL(-1.0, m_one.v());

  const auto d_120{Degrees(120.0)};
  const auto d_m120{Degrees(-120.0)};
  BOOST_CHECK_EQUAL(d_120.v(), d_m120.abs().v());

  BOOST_CHECK_EQUAL(30.0, (Degrees(-155.0) - Degrees(175.0)).v());

  BOOST_CHECK_EQUAL(-120.0, (d_120 + d_120).v());
  BOOST_CHECK_EQUAL(120.0, (d_m120 + d_m120).v());
  BOOST_CHECK_EQUAL(120.0, (d_m120 - d_120).v());

  BOOST_CHECK_EQUAL(-60.0, d_120.opposite().v());
  BOOST_CHECK_EQUAL(60.0, d_m120.opposite().v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_radians) {
  const auto one{Radians(1.0)};
  const auto two{Radians(2.0)};
  const auto m_two{-two};

  BOOST_CHECK_EQUAL(-2.0, m_two.v());
  BOOST_CHECK_EQUAL(two.v(), m_two.abs().v());

  const auto m_one{one + m_two};
  BOOST_CHECK_EQUAL(-1.0, m_one.v());

  const auto result_1 = m_two - two;
  BOOST_CHECK_EQUAL(trig::TAU<double> - 4.0, result_1.v());
  BOOST_CHECK_EQUAL(trig::PI<double> - 4.0, result_1.opposite().v());

  const auto result_2 = two - m_two;
  BOOST_CHECK_EQUAL(4.0 - trig::TAU<double>, result_2.v());
  BOOST_CHECK_EQUAL(4.0 - trig::PI<double>, result_2.opposite().v());

  const auto small_negative{Radians(-EPSILON)};
  BOOST_CHECK_EQUAL(0.0, small_negative.clamp(one).v());
  const auto zero{Radians(0.0)};
  BOOST_CHECK_EQUAL(0.0, zero.clamp(one).v());
  BOOST_CHECK_EQUAL(1.0, one.clamp(one).v());
  const auto one_epsilon{Radians(1.0 + EPSILON)};
  BOOST_CHECK_EQUAL(1.0, one_epsilon.clamp(one).v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Angle_default_constructor) {
  // Test default constructed values
  const Angle<double> zero;
  BOOST_CHECK(zero.is_valid());
  BOOST_CHECK(zero);
  BOOST_CHECK_EQUAL(0.0, zero.sin().v());
  BOOST_CHECK_EQUAL(1.0, zero.cos().v());
  BOOST_CHECK_EQUAL(0.0, zero.tan().value());
  BOOST_CHECK_EQUAL(0.0, zero.to_degrees().v());
  BOOST_CHECK_EQUAL(0.0, zero.to_radians().v());
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
