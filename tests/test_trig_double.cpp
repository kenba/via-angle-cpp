//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024-2025 Ken Barker
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
/// @file test_trig_double.hpp
/// @brief Contains the via::trig namespace.
//////////////////////////////////////////////////////////////////////////////
#include "via/angle/trig.hpp"
#include <boost/test/unit_test.hpp>
#include <sstream>

using namespace via::trig;

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
constexpr auto CALCULATION_TOLERANCE{101 * EPSILON};
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_trig)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_pi_functions) {
  // test double precision
  BOOST_CHECK_EQUAL(std::acos(-1), std::numbers::pi_v<double>);
  BOOST_CHECK_EQUAL(std::numbers::pi_v<double> / 180.0, deg2rad(1.0));
  BOOST_CHECK_EQUAL(180.0 / std::numbers::pi_v<double>, rad2deg(1.0));

  BOOST_CHECK_CLOSE(std::acos(-1.0L), PI<long double>,
                    4 * std::numeric_limits<long double>::epsilon());
  BOOST_CHECK_EQUAL(std::numbers::pi_v<long double>, PI<long double>);
  BOOST_CHECK_EQUAL(std::numbers::pi_v<long double> / 180.0L, deg2rad(1.0L));
  BOOST_CHECK_EQUAL(180.0L / std::numbers::pi_v<long double>, rad2deg(1.0L));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_UnitNegRange_traits) {
  const auto zero{UnitNegRange<double>()};
  BOOST_CHECK_EQUAL(0.0, zero.v());
  const auto one{UnitNegRange<double>(1)};
  const auto minus_one{UnitNegRange<double>(-1)};

  BOOST_CHECK_EQUAL(one, one);
  BOOST_CHECK(minus_one < one);
  BOOST_CHECK(minus_one <= one);

  BOOST_CHECK(minus_one != one);
  BOOST_CHECK(one > minus_one);
  BOOST_CHECK(one >= minus_one);

  BOOST_CHECK_EQUAL(one, minus_one.abs());

  std::ostringstream os;
  os << one;
  BOOST_CHECK(std::string("1") == os.str());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_UnitNegRange_clamp) {
  BOOST_CHECK_EQUAL(-1.0, UnitNegRange<double>::clamp(-1.0 - EPSILON).v());
  BOOST_CHECK_EQUAL(-1.0, UnitNegRange<double>::clamp(-1.0).v());
  BOOST_CHECK_EQUAL(1.0, UnitNegRange<double>::clamp(1.0).v());
  BOOST_CHECK_EQUAL(1.0, UnitNegRange<double>::clamp(1.0 + EPSILON).v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_trig_functions) {
  const UnitNegRange cos_60{0.5};
  const auto sin_60{swap_sin_cos(cos_60)};

  const auto sin_120{sin_60};
  const auto cos_120{cosine_from_sine(sin_120, -1.0)};

  const auto zero{cosine_from_sine(UnitNegRange(1.0), -1.0)};
  BOOST_CHECK_EQUAL(0.0, zero.v());
  BOOST_CHECK(!std::signbit(zero.v()));

  const auto recip_sq_epsilon{1.0 / via::trig::SQ_EPSILON<double>};

  const auto sin_msq_epsilon{UnitNegRange(-via::trig::SQ_EPSILON<double>)};
  BOOST_CHECK_EQUAL(-recip_sq_epsilon, csc(sin_msq_epsilon).value());
  BOOST_CHECK_EQUAL(-recip_sq_epsilon, sec(sin_msq_epsilon).value());

  const auto cos_msq_epsilon{swap_sin_cos(sin_msq_epsilon)};
  BOOST_CHECK_EQUAL(1.0, sec(cos_msq_epsilon).value());
  BOOST_CHECK_EQUAL(1.0, csc(cos_msq_epsilon).value());

  BOOST_CHECK_EQUAL(-via::trig::SQ_EPSILON<double>,
                    tan(sin_msq_epsilon, cos_msq_epsilon).value());
  BOOST_CHECK_EQUAL(-recip_sq_epsilon,
                    cot(sin_msq_epsilon, cos_msq_epsilon).value());

  // Test differences
  BOOST_CHECK_CLOSE(sin_60.v(), sine_diff(sin_120, cos_120, sin_60, cos_60).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(cos_60.v(),
                    cosine_diff(sin_120, cos_120, sin_60, cos_60).v(),
                    2 * CALCULATION_TOLERANCE);

  // Test sums
  BOOST_CHECK_CLOSE(sin_120.v(), sine_sum(sin_60, cos_60, sin_60, cos_60).v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(cos_120.v(), cosine_sum(sin_60, cos_60, sin_60, cos_60).v(),
                    2 * CALCULATION_TOLERANCE);

  BOOST_CHECK_CLOSE(sin_60.v(), std::sqrt(sq_sine_half(cos_120)),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(cos_60.v(), std::sqrt(sq_cosine_half(cos_120)),
                    CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_radians_conversion) {
  // -π/6 radians round trip
  const auto [sin_m30, cos_m30]{sincos(-PI_6<double>)};
  BOOST_CHECK_EQUAL(-0.5, sin_m30.v());
  BOOST_CHECK_EQUAL(COS_30_DEGREES<double>, cos_m30.v());
  BOOST_CHECK_EQUAL(-PI_6<double>, arctan2(sin_m30, cos_m30));

  // π/3 radians round trip
  const auto [sin_60, cos_60]{sincos(PI_3<double>)};
  // Not exactly correct because PI is an irrational number
  // BOOST_CHECK_EQUAL(COS_30_DEGREES<double>, sin_60.v());
  BOOST_CHECK_CLOSE(COS_30_DEGREES<double>, sin_60.v(),
                    2 * CALCULATION_TOLERANCE);
  // BOOST_CHECK_EQUAL(0.5, cos_60.v());
  BOOST_CHECK_CLOSE(0.5, cos_60.v(), 2 * CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(PI_3<double>, arctan2(sin_60, cos_60));

  // -π radians round trip to +π radians
  const auto [sin_180, cos_180]{sincos(-PI<double>)};
  BOOST_CHECK_EQUAL(0.0, sin_180.v());
  BOOST_CHECK_EQUAL(-1.0, cos_180.v());
  BOOST_CHECK_EQUAL(PI<double>, arctan2(sin_180, cos_180));

  // π - π/4 radians round trip
  const auto [sin_135, cos_135]{sincos_diff(PI<double>, PI_4<double>)};
  BOOST_CHECK_EQUAL(SQRT1_2<double>, sin_135.v());
  BOOST_CHECK_EQUAL(-SQRT1_2<double>, cos_135.v());
  BOOST_CHECK_EQUAL(PI<double> - PI_4<double>, arctan2(sin_135, cos_135));

  // 6*π - π/3 radians round trip
  const auto [sin_m60, cos_m60]{sincos_diff(3 * TAU<double>, PI_3<double>)};
  // Not exactly correct because PI is an irrational number
  // BOOST_CHECK_EQUAL(-COS_30_DEGREES<double>, sin_m60.v());
  BOOST_CHECK_CLOSE(-COS_30_DEGREES<double>, sin_m60.v(),
                    2 * CALCULATION_TOLERANCE);
  // BOOST_CHECK_EQUAL(0.5, cos_m60.v());
  BOOST_CHECK_CLOSE(0.5, cos_m60.v(), 2 * CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(-PI_3<double>, arctan2(sin_m60, cos_m60));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_degrees_conversion) {
  // -30 degrees round trip
  const auto [sin_m30, cos_m30]{sincosd<double>(-30)};
  BOOST_CHECK_EQUAL(-0.5, sin_m30.v());
  BOOST_CHECK_EQUAL(COS_30_DEGREES<double>, cos_m30.v());
  BOOST_CHECK_EQUAL(-30.0, arctan2d(sin_m30, cos_m30));

  // 60 degrees round trip
  const auto [sin_60, cos_60]{sincosd<double>(60)};
  BOOST_CHECK_EQUAL(COS_30_DEGREES<double>, sin_60.v());
  BOOST_CHECK_EQUAL(0.5, cos_60.v());
  BOOST_CHECK_EQUAL(60.0, arctan2d(sin_60, cos_60));

  // -180 degrees round trip to +180 degrees
  const auto [sin_180, cos_180]{sincosd<double>(-180)};
  BOOST_CHECK_EQUAL(0.0, sin_180.v());
  BOOST_CHECK_EQUAL(-1.0, cos_180.v());
  BOOST_CHECK_EQUAL(180.0, arctan2d(sin_180, cos_180));

  // 180 - 45 degrees round trip
  const auto [sin_135, cos_135]{sincosd_diff<double>(180, 45)};
  BOOST_CHECK_EQUAL(SQRT1_2<double>, sin_135.v());
  BOOST_CHECK_EQUAL(-SQRT1_2<double>, cos_135.v());
  BOOST_CHECK_EQUAL(135.0, arctan2d(sin_135, cos_135));

  // 1080 - 60 degrees round trip
  const auto [sin_m60, cos_m60]{sincosd_diff<double>(1080, 60)};
  BOOST_CHECK_EQUAL(-COS_30_DEGREES<double>, sin_m60.v());
  BOOST_CHECK_EQUAL(0.5, cos_m60.v());
  BOOST_CHECK_EQUAL(-60.0, arctan2d(sin_m60, cos_m60));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_adjacent_length) {
  // length == hypotenuse
  BOOST_CHECK_EQUAL(0.0, calculate_adjacent_length(5.0, 5.0));

  // length == 0.0
  BOOST_CHECK_EQUAL(5.0, calculate_adjacent_length(0.0, 5.0));

  // length > hypotenuse
  BOOST_CHECK_EQUAL(0.0, calculate_adjacent_length(6.0, 5.0));

  // 3, 4, 5 triangle
  BOOST_CHECK_EQUAL(3.0, calculate_adjacent_length(4.0, 5.0));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_spherical_adjacent_length) {
  // length == hypotenuse
  BOOST_CHECK_EQUAL(0.0, spherical_adjacent_length(deg2rad(5.0), deg2rad(5.0)));

  // length == 0
  BOOST_CHECK_EQUAL(deg2rad(5.0), spherical_adjacent_length(0.0, deg2rad(5.0)));

  // length > hypotenuse
  BOOST_CHECK_EQUAL(0.0, spherical_adjacent_length(deg2rad(6.0), deg2rad(5.0)));

  // 3, 4, 5 triangle
  BOOST_CHECK_CLOSE(deg2rad(3.0),
                    spherical_adjacent_length(deg2rad(4.0), deg2rad(5.0)), 0.1);

  // small 3, 4, 5 triangle
  BOOST_CHECK_CLOSE(
      deg2rad(3.0) / 60,
      spherical_adjacent_length(deg2rad(4.0) / 60, deg2rad(5.0) / 60), 2.26e-5);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_spherical_hypotenuse_length) {
  constexpr double zero{0.0};
#ifdef _MSC_VER
  const double three{deg2rad(3.0)};
  const double four{deg2rad(4.0)};
  const double five{deg2rad(5.0)};
#else
  constexpr double three{deg2rad(3.0)};
  constexpr double four{deg2rad(4.0)};
  constexpr double five{deg2rad(5.0)};
#endif

  // Check abort if either a or b are negative
  // Negative length a
  // BOOST_CHECK_EQUAL(three, spherical_hypotenuse_length(-four, three));
  // Negative length b
  // BOOST_CHECK_EQUAL(four, spherical_hypotenuse_length(four, -three));

  // Zero length a
  BOOST_CHECK_EQUAL(three, spherical_hypotenuse_length(zero, three));
  // Zero length b
  BOOST_CHECK_EQUAL(four, spherical_hypotenuse_length(four, zero));
  // Zero length a & b
  BOOST_CHECK_EQUAL(zero, spherical_hypotenuse_length(zero, zero));

  // 3, 4, 5 triangles, note 5 degrees is 0.08726646259971647 radians
  constexpr double result{0.087240926337265545};
  BOOST_CHECK_EQUAL(result, spherical_hypotenuse_length(four, three));
  BOOST_CHECK_EQUAL(result, spherical_hypotenuse_length(three, four));

  // small 3, 4, 5 triangles
  BOOST_CHECK_CLOSE(five / 60,
                    spherical_hypotenuse_length(four / 60, three / 60), 1.0e-5);
  BOOST_CHECK_CLOSE(five / 60,
                    spherical_hypotenuse_length(three / 60, four / 60), 1.0e-5);
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
