#pragma once

//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024 Ken Barker. All Rights Reserved.
// (ken dot barker at via-technology dot co dot uk)
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//////////////////////////////////////////////////////////////////////////////
/// @file two_sum.hpp
/// @brief Contains the two_sum function.
//////////////////////////////////////////////////////////////////////////////
#include <concepts>
#include <tuple>

/// A function to mark a variable NOT to be optimised by the compiler.
/// See: [Enforcing statement order in
/// C++](https://stackoverflow.com/a/38025837)
template <typename T>
#ifdef _MSC_VER
void do_not_optimize(const T& value) {
}
#else
__attribute__((always_inline)) inline void do_not_optimize(const T &value) {
  asm volatile("" : "+m"(const_cast<T &>(value)));
}
#endif

namespace via {
/// Calculate floating-point sum and error.
/// The [2Sum](https://en.wikipedia.org/wiki/2Sum) algorithm.
///
/// @param a, b the floating-point numbers to add.
///
/// returns (a + b) and the floating-point error: $t = a + b - (a \oplus b)$
/// so: $a+b=s+t$.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto two_sum(const T a, const T b) noexcept -> std::tuple<T, T> {
  const T s{a + b};
  const T a_prime{s - b};
  do_not_optimize(&a_prime);
  const T b_prime{s - a_prime};
  do_not_optimize(&b_prime);
  const T delta_a{a - a_prime};
  do_not_optimize(&delta_a);
  const T delta_b{b - b_prime};
  do_not_optimize(&delta_b);
  return {s, delta_a + delta_b};
}
} // namespace via
