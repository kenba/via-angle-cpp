//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2026 Ken Barker
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
/// @file vector2.hpp
/// @brief Contains 2D vector functions.
//////////////////////////////////////////////////////////////////////////////
// #include <cmath>

namespace via {
namespace vector2 {

/// 2D vector double dot product function: a . b.
///
/// @param a_0, a_1 the first vector values.
/// @param b_0, b_1 the second vector values.
///
/// @return the dot product of the 2D vectors.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto dot_product(const T a_0, const T a_1, const T b_0,
                           const T b_1) noexcept -> T {
  // return std::fma(a_1, b_1, a_0 * b_0);
  return a_0 * b_0 + a_1 * b_1;
}

/// 2D vector double perp product function: a x b.
///
/// @param a_0, a_1 the first vector values.
/// @param b_0, b_1 the second vector values.
///
/// @return the perp product of the 2D vectors.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto perp_product(const T a_0, const T a_1, const T b_0,
                            const T b_1) noexcept -> T {
  return dot_product(a_0, a_1, b_1, -b_0);
}

} // namespace vector2
} // namespace via
