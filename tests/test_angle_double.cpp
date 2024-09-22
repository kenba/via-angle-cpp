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
#include "via/trig.hpp"
#include <boost/test/unit_test.hpp>

using namespace via;

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
constexpr auto CALCULATION_TOLERANCE{101 * EPSILON};
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_angle)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_degrees) {
  const auto one{Degrees(1.0)};
  const auto two{Degrees(2.0)};
  const auto m_two{-two};

  BOOST_CHECK_EQUAL(-2.0, m_two.v());
  BOOST_CHECK_EQUAL(two, m_two.abs());

  const auto m_one{one + m_two};
  BOOST_CHECK_EQUAL(-1.0, m_one.v());

  const auto d_120{Degrees(120.0)};
  const auto d_m120{Degrees(-120.0)};
  BOOST_CHECK_EQUAL(d_120, d_m120.abs());

  BOOST_CHECK_EQUAL(30.0, (Degrees(-155.0) - Degrees(175.0)).v());

  BOOST_CHECK_EQUAL(d_m120, (d_120 + d_120));
  BOOST_CHECK_EQUAL(d_120, (d_m120 + d_m120));
  BOOST_CHECK_EQUAL(d_120, (d_m120 - d_120));

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
  BOOST_CHECK_EQUAL(two, m_two.abs());

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
  BOOST_CHECK(!zero.csc());
  BOOST_CHECK_EQUAL(1.0, zero.sec().value());
  BOOST_CHECK(!zero.cot());
  BOOST_CHECK_EQUAL(0.0, zero.to_degrees().v());
  BOOST_CHECK_EQUAL(0.0, zero.to_radians().v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Angle_y_x_constructor) {
  const Angle<double> one(1.0, 0.0);
  BOOST_CHECK(one.is_valid());
  BOOST_CHECK(one);
  BOOST_CHECK_EQUAL(1.0, one.sin().v());
  BOOST_CHECK_EQUAL(0.0, one.cos().v());
  BOOST_CHECK(!one.tan());
  BOOST_CHECK_EQUAL(1.0, one.csc().value());
  BOOST_CHECK(!one.sec());
  BOOST_CHECK_EQUAL(0.0, one.cot().value());
  BOOST_CHECK_EQUAL(90.0, one.to_degrees().v());
  BOOST_CHECK_EQUAL(trig::PI_2<double>, one.to_radians().v());

  const Angle angle_m45(-EPSILON, EPSILON);
  BOOST_CHECK_CLOSE(-trig::SQRT1_2<double>, angle_m45.sin().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(trig::SQRT1_2<double>, angle_m45.cos().v(),
                    CALCULATION_TOLERANCE);

  BOOST_CHECK(angle_m45 < one);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Angle_conversion) {
  const Angle<double> zero;

  const Angle too_small(-EPSILON / 2, EPSILON / 2);
  BOOST_CHECK(too_small.is_valid());
  BOOST_CHECK_EQUAL(zero, too_small);

  const Angle small(Radians(-3.35e7 * EPSILON));
  BOOST_CHECK(small.is_valid());
  BOOST_CHECK_EQUAL(-3.35e7 * EPSILON, small.sin().v());
  BOOST_CHECK_EQUAL(1.0, small.cos().v());

  const Angle angle_30(Radians(trig::PI_3<double>),
                       Radians(trig::PI_6<double>));
  BOOST_CHECK(angle_30.is_valid());
  BOOST_CHECK_EQUAL(0.5, angle_30.sin().v());
  BOOST_CHECK_EQUAL(trig::SQRT3<double> / 2, angle_30.cos().v());
  BOOST_CHECK_EQUAL(Degrees(30.0), angle_30.to_degrees());
  BOOST_CHECK_EQUAL(Radians(trig::PI_6<double>), angle_30.to_radians());

  const Angle angle_45(Radians(trig::PI_4<double>));
  BOOST_CHECK(angle_45.is_valid());
  BOOST_CHECK_EQUAL(trig::SQRT1_2<double>, angle_45.sin().v());
  BOOST_CHECK_EQUAL(trig::SQRT1_2<double>, angle_45.cos().v());
  BOOST_CHECK_EQUAL(Degrees(45.0), angle_45.to_degrees());
  BOOST_CHECK_EQUAL(Radians(trig::PI_4<double>), angle_45.to_radians());

  const Angle angle_m45(Degrees(-45.0));
  BOOST_CHECK(angle_m45.is_valid());
  BOOST_CHECK_EQUAL(-trig::SQRT1_2<double>, angle_m45.sin().v());
  BOOST_CHECK_EQUAL(trig::SQRT1_2<double>, angle_m45.cos().v());
  BOOST_CHECK_EQUAL(Degrees(-45.0), angle_m45.to_degrees());
  BOOST_CHECK_EQUAL(Radians(-trig::PI_4<double>), angle_m45.to_radians());

  const Angle angle_60(Degrees(-140.0), Degrees(160.0));
  BOOST_CHECK(angle_60.is_valid());
  BOOST_CHECK_EQUAL(trig::SQRT3<double> / 2, angle_60.sin().v());
  BOOST_CHECK_EQUAL(0.5, angle_60.cos().v());
  BOOST_CHECK_EQUAL(Degrees(60.0), angle_60.to_degrees());
  // Fails because PI is irrational
  // BOOST_CHECK_EQUAL(Radians(trig::PI_3<double>), angle_60.to_radians());
  BOOST_CHECK_CLOSE(trig::PI_3<double>, angle_60.to_radians().v(),
                    CALCULATION_TOLERANCE);

  const Angle angle_d30(Degrees(-155.0), Degrees(175.0));
  BOOST_CHECK(angle_d30.is_valid());
  BOOST_CHECK_EQUAL(0.5, angle_d30.sin().v());
  BOOST_CHECK_EQUAL(trig::SQRT3<double> / 2, angle_d30.cos().v());
  BOOST_CHECK_EQUAL(Degrees(30.0), angle_d30.to_degrees());
  BOOST_CHECK_EQUAL(Radians(trig::PI_6<double>), angle_d30.to_radians());

  const Angle angle_120(Degrees(120.0));
  BOOST_CHECK(angle_120.is_valid());
  BOOST_CHECK_EQUAL(trig::SQRT3<double> / 2, angle_120.sin().v());
  BOOST_CHECK_EQUAL(-0.5, angle_120.cos().v());
  BOOST_CHECK_EQUAL(Degrees(120.0), angle_120.to_degrees());
  BOOST_CHECK_EQUAL(Radians(2 * trig::PI_3<double>), angle_120.to_radians());

  const Angle angle_m120(Degrees(-120.0));
  BOOST_CHECK(angle_m120.is_valid());
  BOOST_CHECK_EQUAL(-trig::SQRT3<double> / 2, angle_m120.sin().v());
  BOOST_CHECK_EQUAL(-0.5, angle_m120.cos().v());
  BOOST_CHECK_EQUAL(Degrees(-120.0), angle_m120.to_degrees());
  BOOST_CHECK_EQUAL(Radians(-2 * trig::PI_3<double>), angle_m120.to_radians());

  const Angle angle_m140(Degrees(-140.0));
  BOOST_CHECK(angle_m140.is_valid());
  BOOST_CHECK_CLOSE(-0.6427876096865393, angle_m140.sin().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(-0.7660444431189781, angle_m140.cos().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(Degrees(-140.0), angle_m140.to_degrees());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_angle_maths) {
  const Angle degrees_30(Degrees(30.0));
  const Angle degrees_60(Degrees(60.0));
  const Angle degrees_120(Degrees(120.0));
  const Angle degrees_m120(Degrees(-120.0));

  BOOST_CHECK(degrees_120 < degrees_m120);
  BOOST_CHECK_EQUAL(degrees_120, degrees_m120.abs());
  BOOST_CHECK_EQUAL(degrees_60, degrees_m120.opposite());
  BOOST_CHECK_EQUAL(degrees_120, degrees_30.quarter_turn_cw());
  BOOST_CHECK_EQUAL(degrees_30, degrees_120.quarter_turn_ccw());
  BOOST_CHECK_EQUAL(degrees_60, degrees_120.negate_cos());

  BOOST_CHECK_EQUAL(Degrees(120.0), (degrees_m120 - degrees_120).to_degrees());
  BOOST_CHECK_EQUAL(Degrees(-120.0), (degrees_120 + degrees_120).to_degrees());
  BOOST_CHECK_EQUAL(Degrees(120.0), degrees_60.x2().to_degrees());
  BOOST_CHECK_EQUAL(Degrees(-120.0), degrees_120.x2().to_degrees());
  BOOST_CHECK_EQUAL(Degrees(-60.0), degrees_m120.half().to_degrees());
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
