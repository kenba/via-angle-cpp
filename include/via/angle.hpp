#pragma once

//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024 Ken Barker
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
//////////////////////////////////////////////////////////////////////////////
#include "angle/trig.hpp"

namespace via {
/// The difference of two angles in degrees, reduced to +/-180 degrees.
/// This is a refactor of GeographicLib::AngDiff.
/// @param x, y the angles in degrees
/// @return the difference and error of the angles.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto degrees_diff(const T x, const T y) noexcept -> std::tuple<T, T> {
  auto [a, b] = two_sum(std::remainder(-x, T(360)), std::remainder(y, T(360)));
  auto [d, e] = two_sum(std::remainder(a, T(360)), b);

  // Fix the sign if d = -180, 0, 180.
  if (d == 0 || std::abs(d) == T(360))
    d = std::copysign(d, (e == 0) ? y - x : -e);

  return (d, e);
}

/// The Degrees type.
template <typename T>
  requires std::floating_point<T>
class Degrees {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  T v_;
#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif
  /// Constructor
  constexpr explicit Degrees(const T value) noexcept : v_{value} {}

  /// The accessor for v.
  [[nodiscard("Pure Function")]]
  constexpr auto v() const noexcept -> T {
    return v_;
  }

  /// The absolute value of the `Degrees`
  [[nodiscard("Pure Function")]]
  constexpr auto abs() const noexcept -> Degrees<T> {
    return Degrees(std::abs(v_));
  }

  /// The opposite angle on the circle, i.e. +/- 180 degrees.
  [[nodiscard("Pure Function")]]
  constexpr auto opposite() const noexcept -> Degrees<T> {
    return Degrees((0 < v_) ? v_ - 180 : v_ + 180);
  }

  /// The unary minus operator
  [[nodiscard("Pure Function")]]
  constexpr auto operator-() const noexcept -> Degrees<T> {
    return Degrees(T() - v_);
  }

  /// A Python representation of a Degrees type.
  /// I.e.: Degrees(v)
  /// @return a string in Python repr format.
  std::string python_repr() const {
    return "Degrees( " + std::to_string(v_) + " )";
  }
};

/// Degrees equality operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator==(const Degrees<T> &lhs, const Degrees<T> &rhs) noexcept
    -> bool {
  return lhs.v() == rhs.v();
}

/// Degrees ostream << operator
template <typename T>
  requires std::floating_point<T>
constexpr auto operator<<(std::ostream &os, const Degrees<T> &a)
    -> std::ostream & {
  return os << a.v();
}

/// Degrees addition operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator+(const Degrees<T> &lhs, const Degrees<T> &rhs) noexcept
    -> Degrees<T> {
  constexpr T T180{180};
  constexpr T T360{360};
  const auto [s, t]{two_sum(lhs.v(), rhs.v())};
  return Degrees((s <= static_cast<T>(-T180)) ? s + T360 + t
                 : (s > T180)                 ? s - T360 + t
                                              : s);
}

/// Degrees subtraction operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr Degrees<T> operator-(const Degrees<T> &lhs,
                               const Degrees<T> &rhs) noexcept {
  return lhs + -rhs;
}

/// The Radians type.
template <typename T>
  requires std::floating_point<T>
class Radians {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  T v_;
#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif
  /// Constructor
  constexpr explicit Radians(const T value) noexcept : v_{value} {}

  /// Degrees Constructor
  constexpr explicit Radians(const Degrees<T> value) noexcept
      : v_{trig::deg2rad(value.v())} {}

  /// The accessor for v.
  [[nodiscard("Pure Function")]]
  constexpr auto v() const noexcept -> T {
    return v_;
  }

  /// The absolute value of the `Radians`
  [[nodiscard("Pure Function")]]
  constexpr auto abs() const noexcept -> Radians<T> {
    return Radians(std::abs(v_));
  }

  /// The opposite angle on the circle, i.e. +/- PI radians.
  [[nodiscard("Pure Function")]]
  constexpr auto opposite() const noexcept -> Radians<T> {
    return Radians((0 < v_) ? v_ - trig::PI<T> : v_ + trig::PI<T>);
  }

  /// Clamp value into the range: `0 <= v_ <= max_value`.
  [[nodiscard("Pure Function")]]
  constexpr auto clamp(const Radians<T> max_value) const noexcept
      -> Radians<T> {
    return Radians(std::clamp<T>(v_, 0, max_value.v()));
  }

  /// The unary minus operator
  [[nodiscard("Pure Function")]]
  constexpr auto operator-() const noexcept -> Radians<T> {
    return Radians(T() - v_);
  }

  /// A Python representation of a Radians type.
  /// I.e.: Radians(v)
  /// @return a string in Python repr format.
  std::string python_repr() const {
    return "Radians( " + std::to_string(v_) + " )";
  }
};

/// Radians equality operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator==(const Radians<T> &lhs, const Radians<T> &rhs) noexcept
    -> bool {
  return lhs.v() == rhs.v();
}

/// Radians ostream << operator
template <typename T>
  requires std::floating_point<T>
constexpr auto operator<<(std::ostream &os, const Radians<T> &a)
    -> std::ostream & {
  return os << a.v();
}

/// Radians addition operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator+(const Radians<T> &lhs, const Radians<T> &rhs) noexcept
    -> Radians<T> {
  const auto [s, t]{two_sum(lhs.v(), rhs.v())};
  return Radians((s <= -trig::PI<T>) ? s + trig::TAU<T> + t
                 : (s > trig::PI<T>) ? s - trig::TAU<T> + t
                                     : s);
}

/// Radians subtraction operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator-(const Radians<T> &lhs, const Radians<T> &rhs) noexcept
    -> Radians<T> {
  return lhs + -rhs;
}

/// The Angle represents an angle by it's sine and cosine components.
/// @invariant sin() * sin() + cos() * cos() = 1
template <typename T>
  requires std::floating_point<T>
class Angle {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  trig::UnitNegRange<T> sin_ =
      trig::UnitNegRange<T>(0); ///< The sine of the angle.
  trig::UnitNegRange<T> cos_ =
      trig::UnitNegRange<T>(1); ///< The cosine of the angle.

#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif

  /// Default constructor. Zero angle, i.e. North.
  constexpr Angle() noexcept = default;

  /// Element constructor. Sets the values.
  /// @pre s * s + c * c = 1
  /// @post is_valid() == true.
  /// @param s the sine of the angle.
  /// @param c the cosine of the angle.
  constexpr Angle(const trig::UnitNegRange<T> s, const trig::UnitNegRange<T> c)
      : sin_{s}, cos_{c} {
#ifndef PYBIND11_VERSION_MAJOR
    Ensures(is_valid());
#endif
  }

  /// Tuple Constructor
  constexpr Angle(
      const std::tuple<trig::UnitNegRange<T>, trig::UnitNegRange<T>> &t)
      : Angle(std::get<0>(t), std::get<1>(t)) {}

  /// Degrees Constructor
  constexpr explicit Angle(const Degrees<T> degrees)
      : Angle(trig::sincosd(degrees.v())) {}

  /// Degrees difference Constructor
  constexpr explicit Angle(const Degrees<T> a, const Degrees<T> b)
      : Angle(trig::sincosd_diff(a.v(), b.v())) {}

  /// Radians Constructor
  constexpr explicit Angle(const Radians<T> radians)
      : Angle(trig::sincos(radians.v())) {}

  /// Radians difference Constructor
  constexpr explicit Angle(const Radians<T> a, const Radians<T> b)
      : Angle(trig::sincos_diff(a.v(), b.v())) {}

  /// Construct an Angle from y and x values.
  /// Normalizes the values.
  /// @post is_valid() == true.
  /// @param y the y coordinate.
  /// @param x the x coordinate.
  static constexpr auto from_y_x(const T y, const T x) -> Angle<T> {
    const T length{std::hypot(y, x)};

    if (length < std::numeric_limits<T>::epsilon()) {
      return Angle<T>();
    } else {
      return Angle(trig::UnitNegRange<T>::clamp(y / length),
                   trig::UnitNegRange<T>::clamp(x / length));
    }
  }

  /// Function to determine whether the Angle is valid.
  [[nodiscard("Pure Function")]]
  constexpr auto is_valid() const noexcept -> bool {
    constexpr T P_EPSILON{32 * std::numeric_limits<T>::epsilon()};

    const T sq_length{sin_.v() * sin_.v() + cos_.v() * cos_.v()};
    return ((1 - P_EPSILON) <= sq_length) && (sq_length <= (1 + P_EPSILON));
  }

  /// Boolean operator to test whether the Angle is valid.
  [[nodiscard("Pure Function")]]
  explicit constexpr operator bool() const noexcept {
    return is_valid();
  }

  /// Accessor for the sine of the Angle.
  [[nodiscard("Pure Function")]]
  constexpr auto sin() const noexcept -> trig::UnitNegRange<T> {
    return sin_;
  }

  /// Accessor for the cosine of the Angle.
  [[nodiscard("Pure Function")]]
  constexpr auto cos() const noexcept -> trig::UnitNegRange<T> {
    return cos_;
  }

  /// The tangent of the Angle.
  /// returns the tangent if `|cos| >= SQ_EPSILON`
  [[nodiscard("Pure Function")]]
  constexpr auto tan() const noexcept -> std::optional<T> {
    return trig::tan(sin_, cos_);
  }

  /// The cosecant of the Angle.
  /// returns the cosecant if `|sin| >= SQ_EPSILON`
  [[nodiscard("Pure Function")]]
  constexpr auto csc() const noexcept -> std::optional<T> {
    return trig::csc(sin_);
  }

  /// The secant of the Angle.
  /// returns the secant if `|cos| >= SQ_EPSILON`
  [[nodiscard("Pure Function")]]
  constexpr auto sec() const noexcept -> std::optional<T> {
    return trig::sec(cos_);
  }

  /// The cotangent of the Angle.
  /// returns the cotangent if `|sin| >= SQ_EPSILON`
  [[nodiscard("Pure Function")]]
  constexpr auto cot() const noexcept -> std::optional<T> {
    return trig::cot(sin_, cos_);
  }

  /// Convert the Angle to Degrees.
  [[nodiscard("Pure Function")]]
  constexpr auto to_degrees() const noexcept -> Degrees<T> {
    return Degrees<T>(trig::arctan2d(sin_, cos_));
  }

  /// Convert the Angle to Radiuans.
  [[nodiscard("Pure Function")]]
  constexpr auto to_radians() const noexcept -> Radians<T> {
    return Radians<T>(trig::arctan2(sin_, cos_));
  }

  /// The absolute value of the Angle, i.e. the Angle with a positive sine.
  /// @return the absolute value angle, the positive angle.
  [[nodiscard("Pure Function")]]
  constexpr auto abs() const noexcept -> Angle<T> {
    return Angle(sin_.abs(), cos_);
  }

  /// The opposite angle on the circle, i.e. +/- 180 degrees.
  /// @return the opposite angle, i.e. swap both signs.
  [[nodiscard("Pure Function")]]
  constexpr auto opposite() const noexcept -> Angle<T> {
    return Angle(-sin_, -cos_);
  }

  /// A quarter turn clockwise around the circle, i.e. + 90 degrees.
  /// @return the angle rotated 90 degrees clockwise.
  [[nodiscard("Pure Function")]]
  constexpr auto quarter_turn_cw() const noexcept -> Angle<T> {
    return Angle(cos_, -sin_);
  }

  /// A quarter turn counter clockwise around the circle, i.e. - 90 degrees.
  /// @return the angle rotated 90 degrees counter clockwise.
  [[nodiscard("Pure Function")]]
  constexpr auto quarter_turn_ccw() const noexcept -> Angle<T> {
    return Angle(-cos_, sin_);
  }

  /// Negate the cosine of the Angle.
  /// I.e. `PI` - `angle.radians()` for positive angles,
  ///      `angle.radians()` + `PI` for negative angles
  /// @return PI - angle.
  [[nodiscard("Pure Function")]]
  constexpr auto negate_cos() const noexcept -> Angle<T> {
    return Angle(sin_, -cos_);
  }

  /// Double the Angle.
  /// See: [Double-angle
  /// formulae](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Double-angle_formulae)
  /// @return two times the angle.
  [[nodiscard("Pure Function")]]
  constexpr auto x2() const noexcept -> Angle<T> {
    return Angle(trig::UnitNegRange<T>::clamp(2 * sin_.v() * cos_.v()),
                 trig::UnitNegRange<T>::clamp((cos_.v() - sin_.v()) *
                                              (cos_.v() + sin_.v())));
  }

  /// Calculate half the angle.
  /// @return half the angle.
  [[nodiscard("Pure Function")]]
  constexpr auto half() const noexcept -> Angle<T> {
    return Angle(
        trig::UnitNegRange<T>::clamp(
            std::copysign(std::sqrt(trig::sq_sine_half(cos_)), sin_.v())),
        trig::UnitNegRange<T>::clamp(std::sqrt(trig::sq_cosine_half(cos_))));
  }

  /// Unary minus operator.
  /// @return the negative angle, i.e. just swap the sine's sign.
  [[nodiscard("Pure Function")]]
  constexpr auto operator-() const noexcept -> Angle<T> {
    return Angle(-sin_, cos_);
  }

  /// Addition operator.
  [[nodiscard("Pure Function")]]
  constexpr friend auto operator+(const Angle<T> a, const Angle<T> &b) noexcept
      -> Angle<T> {
    return Angle(sine_sum(a, b), cosine_sum(a, b));
  }

  /// Subtraction operator.
  [[nodiscard("Pure Function")]]
  constexpr friend auto operator-(const Angle<T> a, const Angle<T> &b) noexcept
      -> Angle<T> {
    return Angle(sine_diff(a, b), cosine_diff(a, b));
  }

  /// The spaceship operator
  /// It compares whether an `Angle` is clockwise of the other `Angle` on the
  /// unit circle.
  [[nodiscard("Pure Function")]]
  constexpr std::partial_ordering operator<=>(const Angle<T> &other) const {
    const auto delta = sine_diff(*this, other);
    return delta <=> trig::UnitNegRange<T>(0);
  }

  /// A Python representation of an Angle type.
  /// I.e.: Angle([sin, cos])
  /// @return a string in Python repr format.
  std::string python_repr() const {
    return "Angle([ " + std::to_string(sin_.v()) + ", " +
           std::to_string(cos_.v()) + " ])";
  }
};

/// Angle equality operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator==(const Angle<T> &lhs, const Angle<T> &rhs) noexcept
    -> bool {
  return lhs.sin() == rhs.sin() && lhs.cos() == rhs.cos();
}

/// Angle ostream << operator
template <typename T>
  requires std::floating_point<T>
constexpr auto operator<<(std::ostream &os, const Angle<T> &a)
    -> std::ostream & {
  return os << '(' << a.sin() << ',' << a.cos() << ')';
}

/// Calculate the sine of the sum of two angles.
/// @param a, b the angles.
/// @return sin(a + b)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sine_sum(const Angle<T> a, const Angle<T> b) noexcept
    -> trig::UnitNegRange<T> {
  return trig::sine_sum(a.sin(), a.cos(), b.sin(), b.cos());
}

/// Calculate the sine of the difference of two angles.
/// @param a, b the angles.
/// @return sin(a - b)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sine_diff(const Angle<T> a, const Angle<T> b) noexcept
    -> trig::UnitNegRange<T> {
  return trig::sine_diff(a.sin(), a.cos(), b.sin(), b.cos());
}

/// Calculate the cosine of the sum of two angles.
/// @param a, b the angles.
/// @return cos(a + b)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto cosine_sum(const Angle<T> a, const Angle<T> b) noexcept
    -> trig::UnitNegRange<T> {
  return trig::cosine_sum(a.sin(), a.cos(), b.sin(), b.cos());
}

/// Calculate the cosine of the difference of two angles.
/// @param a, b the angles.
/// @return cos(a - b)
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto cosine_diff(const Angle<T> a, const Angle<T> b) noexcept
    -> trig::UnitNegRange<T> {
  return trig::cosine_diff(a.sin(), a.cos(), b.sin(), b.cos());
}
} // namespace via
