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
inline void do_not_optimize(const T &) {
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
#ifdef _MSC_VER
#pragma optimize( "", off )
#endif
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
