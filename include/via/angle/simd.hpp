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
/// @file simd.hpp
/// @brief Contains simd functions.
//////////////////////////////////////////////////////////////////////////////
#include <immintrin.h>

namespace via {
namespace simd {

/// 2D vector double dot product function: a . b.
///
/// @param a, b the __m128d double vectors.
///
/// @return the 2D dot product of the vectors.
[[nodiscard("Pure Function")]]
inline constexpr auto dot2d(const __m128d a, const __m128d b) noexcept
    -> double {
  const auto c{_mm_dp_pd(a, b, 0x33)};
  return _mm_cvtsd_f64(c);
}

/// 2D vector double dot product function: a . b.
///
/// @param a_0, a_1 the first vector values.
/// @param b_0, b_1 the second vector values.
///
/// @return the 2D dot product of the vectors.
[[nodiscard("Pure Function")]]
inline constexpr auto dot_product(const double a_0, const double a_1,
                                  const double b_0, const double b_1) noexcept
    -> double {
  return dot2d(_mm_set_pd(a_1, a_0), _mm_set_pd(b_1, b_0));
}

/// 2D vector double perp product function: a x b.
///
/// @param a_0, a_1 the first vector values.
/// @param b_0, b_1 the second vector values.
///
/// @return the 2D perp product of the vector double values.
[[nodiscard("Pure Function")]]
inline constexpr auto perp_product(const double a_0, const double a_1,
                                   const double b_0, const double b_1) noexcept
    -> double {
  return dot2d(_mm_set_pd(a_1, a_0), _mm_set_pd(-b_0, b_1));
}

} // namespace simd
} // namespace via
