#pragma once

//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024 Ken Barker. All Rights Reserved.
// (ken dot barker at via-technology dot co dot uk)
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//////////////////////////////////////////////////////////////////////////////
/// @file trig.hpp
/// @brief Contains the via::trig namespace.
//////////////////////////////////////////////////////////////////////////////
#include "two_sum.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <optional>
#include <ostream>
#include <string>
#include <type_traits>
#ifndef PYBIND11_VERSION_MAJOR
#include <gsl/assert>
#endif

#ifndef M_PI // M_PI is not part of the C or C++ standards: _USE_MATH_DEFINES
constexpr double M_PI{3.14159265358979323846};
constexpr double M_PI_2{1.57079632679489661923};
#endif

#ifndef M_PIl // M_PIl may be defined on Linux platforms as a GNU Extension.
constexpr long double M_PIl{3.141592653589793238462643383279502884L};
#endif

#ifndef M_PI_3l // from Rust:
                // https://doc.rust-lang.org/stable/std/f64/consts/constant.FRAC_PI_3.html
constexpr long double M_PI_3l{1.04719755119659774615421446109316763L};
#endif

#ifndef M_PI_6l // from Rust:
                // https://doc.rust-lang.org/stable/std/f64/consts/constant.FRAC_PI_6.html
constexpr long double M_PI_6l{0.52359877559829887307710723054658381L};
#endif

#ifndef M_SQRT1_2 // M_SQRT1_2 is not part of the C or C++ standards:
                  // _USE_MATH_DEFINES
constexpr double M_SQRT1_2{0.70710678118654752440};
#endif

#ifndef M_SQRT1_2l // M_SQRT1_2l may be defined on Linux platforms as a GNU
                   // Extension.
constexpr long double M_SQRT1_2l{0.707106781186547524400844362104849039L};
#endif

#ifndef M_SQRT3l // sqrt(3)
constexpr long double M_SQRT3l{1.732050807568877293527446341505872367L};
#endif

namespace via {
/// ε * ε, a very small number.
/// Where: ε, the difference between 1.0 and the next value representable by
/// type T.
template <typename T>
  requires std::floating_point<T>
constexpr T SQ_EPSILON{std::numeric_limits<T>::epsilon() *
                       std::numeric_limits<T>::epsilon()};

namespace trig {
/// PI to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T PI{
    static_cast<T>(std::is_same<T, long double>::value ? M_PIl : M_PI)};

/// 2 * PI to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T TAU{2 * PI<T>};

/// PI / 2 to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T PI_2{PI<T> / 2};

/// PI / 3 to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T PI_3{static_cast<T>(M_PI_3l)};

/// PI / 4 to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T PI_4{PI_2<T> / 2};

/// PI / 6 to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T PI_6{static_cast<T>(M_PI_6l)};

/// 1/√2 to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T SQRT1_2{static_cast<T>(
    std::is_same<T, long double>::value ? M_SQRT1_2l : M_SQRT1_2)};

/// √3 to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T SQRT3{static_cast<T>(M_SQRT3l)};

/// cos(30°) i.e. √3/2 to double or long double precision of depending on T.
template <typename T>
  requires std::floating_point<T>
constexpr T COS_30_DEGREES{static_cast<T>(M_SQRT3l / 2)};

/// Convert a value in degrees to radians.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto deg2rad(const T value) noexcept -> T {
  if (std::abs(value) == static_cast<T>(30))
    return std::copysign(PI_6<T>, value);

  return value * (PI<T> / 180);
}

/// Convert a value in radians to degrees.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto rad2deg(const T value) noexcept -> T {
  return value * (180 / PI<T>);
}

/// The UnitNegRange type.
/// @invariant -1 <= v <= 1
template <typename T>
  requires std::floating_point<T>
class UnitNegRange {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  T v_;
#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif
  /// Constructor
  constexpr explicit UnitNegRange(const T value) noexcept : v_{value} {
#ifndef PYBIND11_VERSION_MAJOR
    Expects((-1 <= v_) && (v_ <= 1));
#endif
  }

  /// Clamp value into the valid range: -1.0 to +1.0 inclusive
  [[nodiscard("Pure Function")]]
  static constexpr UnitNegRange<T> clamp(const T value) noexcept {
    return UnitNegRange(std::clamp<T>(value, -1, 1));
  }

  /// The accessor for v.
  [[nodiscard("Pure Function")]]
  constexpr auto v() const noexcept -> T {
    return v_;
  }

  /// The absolute value of the `UnitNegRange`
  [[nodiscard("Pure Function")]]
  constexpr auto abs() const noexcept -> UnitNegRange<T> {
    return UnitNegRange(std::abs(v_));
  }

  /// The unary minus operator
  [[nodiscard("Pure Function")]]
  constexpr auto operator-() const noexcept -> UnitNegRange<T> {
    return UnitNegRange(-v_);
  }

  /// The spaceship operator
  constexpr std::partial_ordering
  operator<=>(const UnitNegRange<T> &other) const {
    return v_ <=> other.v_;
  }

  /// A Python representation of a UnitNegRange.
  /// I.e.: UnitNegRange(v)
  /// @return a string in Python repr format.
  std::string python_repr() const {
    return "UnitNegRange(" + std::to_string(v_) + ")";
  }

}; // UnitNegRange

/// UnitNegRange equality operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator==(const UnitNegRange<T> &lhs,
                          const UnitNegRange<T> &rhs) noexcept -> bool {
  return lhs.v() == rhs.v();
}

/// UnitNegRange ostream << operator
template <typename T>
  requires std::floating_point<T>
constexpr auto operator<<(std::ostream &os,
                          const UnitNegRange<T> &a) -> std::ostream & {
  return os << a.v();
}

/// Convert a cosine to a sine, or vice versa:
/// return std::sqrt(1 - a) * (1 + a))
/// @pre -1 <= a <= 1
/// @post 0 <= result <= 1
/// @param a the cosine (or sine) of the angle.
/// @return the positive sine (or cosine) of the angle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
swap_sin_cos(const UnitNegRange<T> a) noexcept -> UnitNegRange<T> {
  return UnitNegRange<T>::clamp(std::sqrt((1 - a.v()) * (1 + a.v())));
};

/// Convert a sine to a cosine, or vice versa:
/// return std::sqrt(1 - a) * (1 + a))
/// If the angle > +/- 90 degrees (i.e. PI_2).
/// @pre -1 <= a <= 1
/// @post -1 <= result <= 1, with sign.
/// @param a the cosine (or sine) of the angle.
/// @param sign the sign of the result.
/// @return the cosine (or sine) of the angle.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto cosine_from_sine(const UnitNegRange<T> a,
                                T sign) noexcept -> UnitNegRange<T> {
  return UnitNegRange<T>(std::copysign(swap_sin_cos(a).v(), sign));
}

/// Calculate the sine of an angle in `radians`.
/// Corrects sin ±π/4 to ±1/√2.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sine(const T radians) noexcept {
  const auto radians_abs{std::abs(radians)};
  if (radians_abs == PI_4<T>)
    return UnitNegRange<T>(std::copysign(SQRT1_2<T>, radians));

  return UnitNegRange<T>(std::sin(radians));
}

/// Calculate the cosine of an angle in `radians` using the sine of the angle.
/// Corrects cos π/4 to 1/√2.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto cosine(const T radians,
                      const UnitNegRange<T> sin) noexcept -> UnitNegRange<T> {
  const auto radians_abs{std::abs(radians)};
  if (radians_abs == PI_4<T>)
    return UnitNegRange<T>(std::copysign(SQRT1_2<T>, PI_2<T> - radians_abs));

  return cosine_from_sine(sin, PI_2<T> - radians_abs);
}

/// Assign `sin` and `cos` to the `remquo` quadrant: `q`:
/// - 0: no conversion
/// - 1: rotate 90° clockwise
/// - 2: opposite quadrant
/// - 3: rotate 90° counter-clockwise
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto assign_sin_cos_to_quadrant(
    const UnitNegRange<T> sin, const UnitNegRange<T> cos,
    const int q) noexcept -> std::tuple<UnitNegRange<T>, UnitNegRange<T>> {
  switch (q & 3) {
  case 1:
    return {cos, -sin}; // quarter_turn_cw
  case 2:
    return {-sin, -cos}; // opposite
  case 3:
    return {-cos, sin}; // quarter_turn_ccw
  default:
    return {sin, cos};
  }
}

/// Calculate the sine and cosine of an angle from a value in `Radians`.
/// Note: calculates the cosine of the angle from its sine and overrides both
/// the sine and cosine for π/4 to their correct values: 1/√2
///
/// @param `radians` the angle in `Radians`
///
/// @return sine and cosine of the angle as `UnitNegRange`s.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sincos(const T radians) noexcept
    -> std::tuple<UnitNegRange<T>, UnitNegRange<T>> {
  int q;
  const auto radians_q{std::remquo(radians, PI_2<T>, &q)};

  const auto sin{sine(radians_q)};
  return assign_sin_cos_to_quadrant(sin, cosine(radians_q, sin), q);
}

/// Calculate the sine and cosine of an angle from the difference of a pair of
/// values in `Radians`.
/// Note: calculates the cosine of the angle from its sine and overrides the
/// sine and cosine for π/4 to their correct values: 1/√2
///
/// @param a, b the angles in `Radians`
///
/// @return sine and cosine of a - b as `UnitNegRange`s.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sincos_diff(const T a, const T b) noexcept
    -> std::tuple<UnitNegRange<T>, UnitNegRange<T>> {
  const auto [delta, t]{two_sum(a, -b)};
  int q;
  const auto radians_q{std::remquo(delta, PI_2<T>, &q) + t};

  const auto sin{sine(radians_q)};
  return assign_sin_cos_to_quadrant(sin, cosine(radians_q, sin), q);
}

/// Accurately calculate an angle in `Radians` from its sine and cosine.
///
/// @param sin, cos the sine and cosine of the angle in `UnitNegRange`s.
///
/// @return the angle in `Radians`.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto arctan2(const UnitNegRange<T> sin,
                       const UnitNegRange<T> cos) noexcept -> T {
  const auto sin_abs{sin.abs()};
  const auto cos_abs{cos.abs()};

  // calculate radians in the range 0 <=radians <= π/2
  // default to π/4, i.e. sin_abs == cos_abs
  auto radians{PI_4<T>};
  if (sin_abs != cos_abs)
    radians = (sin_abs < cos_abs)
                  ? std::atan2(sin_abs.v(), cos_abs.v())
                  : PI_2<T> - std::atan2(cos_abs.v(), sin_abs.v());

  // calculate radians in the range 0 <= radians <= π
  // is cos is negative
  if (std::signbit(cos.v()))
    radians = PI<T> - radians;

  // return radians in the range -π <=radians <= π
  return std::copysign(radians, sin.v());
}

/// Calculate the sine and cosine of an angle from a value in `Degrees`.
/// Note: calculates the cosine of the angle from its sine and overrides the
/// sine and cosine for ±30° and ±45° to their correct values.
///
/// @param degrees the angle in `Degrees`
///
/// @return sine and cosine of the angle as `UnitNegRange`s.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sincosd(const T degrees) noexcept
    -> std::tuple<UnitNegRange<T>, UnitNegRange<T>> {
  int q;
  const auto degrees_q{std::remquo(degrees, static_cast<T>(90), &q)};

  // radians_q is radians in range `-FRAC_PI_4..=FRAC_PI_4`
  const auto radians_q = deg2rad(degrees_q);
  const auto sin{sine(radians_q)};
  return assign_sin_cos_to_quadrant(sin, cosine(radians_q, sin), q);
}

/// Calculate the sine and cosine of an angle from the difference of a pair of
/// values in `Degrees`.
/// Note: calculates the cosine of the angle from its sine and overrides the
/// sine and cosine for ±30° and ±45° to their correct values.
///
/// @param  a, b the angles in `Degrees`
///
/// @return sine and cosine of a - b as `UnitNegRange`s.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sincosd_diff(const T a, const T b) noexcept
    -> std::tuple<UnitNegRange<T>, UnitNegRange<T>> {
  const auto [delta, t] = two_sum(a, -b);
  int q;
  const auto degrees_q{std::remquo(delta, static_cast<T>(90), &q) + t};

  // radians_q is radians in range `-FRAC_PI_4..=FRAC_PI_4`
  const auto radians_q = deg2rad(degrees_q);
  const auto sin{sine(radians_q)};
  return assign_sin_cos_to_quadrant(sin, cosine(radians_q, sin), q);
}

/// Accurately calculate a small an angle in `Degrees` from the its sine and
/// cosine. Converts sin of 0.5 to 30°.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto arctan2_degrees(const UnitNegRange<T> sin_abs,
                               const UnitNegRange<T> cos_abs) noexcept -> T {
  return (sin_abs.v() == 0.5) ? static_cast<T>(30)
                              : rad2deg(std::atan2(sin_abs.v(), cos_abs.v()));
}

/// Accurately calculate an angle in `Radians` from its sine and cosine.
///
/// @param sin, cos the sine and cosine of the angle in `UnitNegRange`s.
///
/// @return the angle in `Radians`.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto arctan2d(const UnitNegRange<T> sin,
                        const UnitNegRange<T> cos) noexcept -> T {
  const auto sin_abs{sin.abs()};
  const auto cos_abs{cos.abs()};

  // calculate degrees in the range 0° <=degrees <= 90°
  // default to 90°, i.e. sin_abs == cos_abs
  T degrees{45};
  if (sin_abs != cos_abs)
    degrees = (sin_abs < cos_abs)
                  ? arctan2_degrees(sin_abs, cos_abs)
                  : static_cast<T>(90) - arctan2_degrees(cos_abs, sin_abs);

  // calculate degrees in the range 0° <= degrees <= 180°
  // is cos is negative
  if (std::signbit(cos.v()))
    degrees = static_cast<T>(180) - degrees;

  // return degrees in the range -180° <= degrees <= 180°
  return std::copysign(degrees, sin.v());
}

/// The cosecant of an angle.
///
/// @param sin the sine of the angle.
///
/// returns the cosecant if `|sin| >= SQ_EPSILON`
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto csc(const UnitNegRange<T> sin) noexcept -> std::optional<T> {
  // protect against divide by zero
  if (sin.abs().v() >= SQ_EPSILON<T>)
    return 1.0 / sin.v();

  return {};
}

/// The secant of an angle.
///
/// @param cos the cosine of the angle.
///
/// returns the secant if `|cos| >= SQ_EPSILON`
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sec(const UnitNegRange<T> cos) noexcept -> std::optional<T> {
  return csc(cos);
}

/// The tangent of an angle.
///
/// @param sin the sine of the angle.
/// @param cos the cosine of the angle.
///
/// returns the tangent if `|cos| >= SQ_EPSILON`
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto tan(const UnitNegRange<T> sin,
                   const UnitNegRange<T> cos) noexcept -> std::optional<T> {
  const auto secant{sec(cos)};
  if (secant.has_value())
    return sin.v() * secant.value();

  return {};
}

/// The cotangent of an angle.
///
/// @param sin the sine of the angle.
/// @param cos the cosine of the angle.
///
/// returns the cotangent if `|sin| >= SQ_EPSILON`
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto cot(const UnitNegRange<T> sin,
                   const UnitNegRange<T> cos) noexcept -> std::optional<T> {
  const auto cosecant{csc(sin)};
  if (cosecant.has_value())
    return cos.v() * cosecant.value();

  return {};
}

/// Calculate the sine of the difference of two angles.
/// It is effectively the perp_product of the two angles as vectors.
/// @pre |params| <= 1
/// @post |results| <= 1
/// @param sin_a, cos_a the sine and cosine of angle a.
/// @param sin_b, cos_b the sine and cosine of angle b.
/// @return sin(a - b)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
sine_diff(const UnitNegRange<T> sin_a, const UnitNegRange<T> cos_a,
          const UnitNegRange<T> sin_b,
          const UnitNegRange<T> cos_b) noexcept -> UnitNegRange<T> {
  return UnitNegRange<T>::clamp(sin_a.v() * cos_b.v() - sin_b.v() * cos_a.v());
}

/// Calculate the cosine of the difference of two angles.
/// It is effectively the dot product of the two angles as vectors.
/// @pre |params| <= 1
/// @post |results| <= 1
/// @param sin_a, cos_a the sine and cosine of angle a.
/// @param sin_b, cos_b the sine and cosine of angle b.
/// @return cos(a - b)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
cosine_diff(const UnitNegRange<T> sin_a, const UnitNegRange<T> cos_a,
            const UnitNegRange<T> sin_b,
            const UnitNegRange<T> cos_b) noexcept -> UnitNegRange<T> {
  return UnitNegRange<T>::clamp(cos_a.v() * cos_b.v() + sin_a.v() * sin_b.v());
}

/// Calculate the sine of the sum of two angles.
/// @pre |params| <= 1
/// @post |results| <= 1
/// @param sin_a, cos_a the sine and cosine of angle a.
/// @param sin_b, cos_b the sine and cosine of angle b.
/// @return sin(a + b)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
sine_sum(const UnitNegRange<T> sin_a, const UnitNegRange<T> cos_a,
         const UnitNegRange<T> sin_b,
         const UnitNegRange<T> cos_b) noexcept -> UnitNegRange<T> {
  return sine_diff(sin_a, cos_a, -sin_b, cos_b);
}

/// Calculate the cosine of the sum of two angles.
/// @pre |params| <= 1
/// @post |results| <= 1
/// @param sin_a, cos_a the sine and cosine of angle a.
/// @param sin_b, cos_b the sine and cosine of angle b.
/// @return cos(a + b)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
cosine_sum(const UnitNegRange<T> sin_a, const UnitNegRange<T> cos_a,
           const UnitNegRange<T> sin_b,
           const UnitNegRange<T> cos_b) noexcept -> UnitNegRange<T> {
  return cosine_diff(sin_a, cos_a, -sin_b, cos_b);
}

/// Square of the sine of half the Angle.
/// See: [Half-angle
/// formulae](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Half-angle_formulae)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sq_sine_half(const UnitNegRange<T> cos_a) noexcept -> T {
  return (1 - cos_a.v()) / 2;
}

/// Square of the cosine of half the Angle.
/// See: [Half-angle
/// formulae](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Half-angle_formulae)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sq_cosine_half(const UnitNegRange<T> cos_a) noexcept -> T {
  return (1 + cos_a.v()) / 2;
}

/// Calculates the length of the other side in a right angled
/// triangle, given the length of one side and the hypotenuse.
/// @pre 0 <= length
/// @pre 0 <= hypotenuse
/// @post result <= hypotenuse
/// @param length the length of a side.
/// @param hypotenuse the length of the hypotenuse
/// @return the length of the other side. zero if length >= hypotenuse
/// or the hypotenuse if length <= 0.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_adjacent_length(const T length,
                               const T hypotenuse) noexcept -> T {
#ifndef PYBIND11_VERSION_MAJOR
  Expects((T() <= length) && (T() <= hypotenuse));
#endif
  const T result{(length <= T()) ? hypotenuse
                                 : ((length >= hypotenuse)
                                        ? T()
                                        : std::sqrt((hypotenuse - length) *
                                                    (hypotenuse + length)))};
#ifndef PYBIND11_VERSION_MAJOR
  Ensures(result <= hypotenuse);
#endif
  return result;
}

/// Calculates the length of the other side in a right angled
/// spherical triangle, given the length of one side and the hypotenuse.
/// @pre 0 <= a
/// @pre 0 <= c
/// @post result <= c
/// @param a the length of a side.
/// @param c the length of the hypotenuse
/// @return the length of the other side. zero if a >= c or c if a <= 0.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto spherical_adjacent_length(const T a, const T c) noexcept -> T {
#ifndef PYBIND11_VERSION_MAJOR
  Expects((T() <= a) && (T() <= c));
#endif
  const T result{(a <= T())
                     ? c
                     : ((a >= c) ? T()
                                 : std::acos(std::clamp<T>(
                                       std::cos(c) / std::cos(a), -1, 1)))};
#ifndef PYBIND11_VERSION_MAJOR
  Ensures(result <= c);
#endif
  return result;
}

/// Calculates the length of the hypotenuse in a right angled
/// spherical triangle, given the length of both sides.
/// @pre 0 <= a
/// @pre 0 <= b
/// @post result >= max(a, b)
/// @param a, b the lengths of the sides.
/// @return the length of the hypotenuse. zero if a <= 0 or b <= 0.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto spherical_hypotenuse_length(const T a, const T b) noexcept -> T {
#ifndef PYBIND11_VERSION_MAJOR
  Expects((T() <= a) && (T() <= b));
#endif
  const T result{(a <= T())   ? ((b <= T()) ? T() : b)
                 : (b <= T()) ? a
                              : std::acos(std::clamp<T>(
                                    std::cos(a) * std::cos(b), -1, 1))};
#ifndef PYBIND11_VERSION_MAJOR
  Ensures(result >= std::max(a, b));
#endif
  return result;
}

/// Calculate the length of the adjacent side of a right angled spherical
/// triangle, given the cosine of the angle and length of the hypotenuse.
/// See: [Spherical law of
/// cosines](https://en.wikipedia.org/wiki/Spherical_law_of_cosines)
/// @param cos_angle the cosine of the adjacent angle.
/// @param length the length of the hypotenuse
///
/// @return the length of the opposite side.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto spherical_cosine_rule(const T cos_angle, const T length) noexcept -> T {
  return std::atan(cos_angle * std::tan(length));
}
} // namespace trig
} // namespace via
