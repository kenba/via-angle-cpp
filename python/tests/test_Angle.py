#!/usr/bin/env python

# Copyright (c) 2024 Ken Barker
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#  @file test_Angle
#  @brief Contains unit tests for the via angle classes and functions.

import pytest
import numpy as np
from numpy.testing import assert_almost_equal
from via_angle import Angle, Degrees, Radians, PI_3, SQRT1_2, SQRT3, TAU

EPSILON = np.finfo(float).eps
PI_6 = PI_3 / 2.0

def test_degrees():
    one = Degrees(1.0)
    two = Degrees(2.0)
    m_two = -two

    assert -2.0 == m_two.v()
    assert two == abs(m_two)

    m_one = one + m_two
    assert -1.0 == m_one.v()
    assert one == abs(m_one)

    d_120 = Degrees(120.0)
    d_m120 = -d_120
    assert d_120 == abs(d_m120)

    assert Degrees(30.0) == Degrees(-155.0) - Degrees(175.0)

    assert d_m120 == d_120 + d_120
    assert d_120 == d_m120 + d_m120
    assert d_120 == d_m120 - d_120

    assert Degrees(-60.0) == d_120.opposite()
    assert Degrees(60.0) == d_m120.opposite()

def test_radians():
    one = Radians(1.0)
    two = Radians(2.0)
    m_two = -two

    assert -2.0 == m_two.v()
    assert two == abs(m_two)

    m_one = one + m_two
    assert -1.0 == m_one.v()
    assert one == abs(m_one)

    result_1 = m_two - two
    assert TAU - 4.0 == result_1.v()
    assert np.pi - 4.0 == result_1.opposite().v()

    result_2 = two - m_two
    assert 4.0 - TAU == result_2.v()
    assert 4.0 - np.pi == result_2.opposite().v()

def test_Angle_default_constructor():
    """Test default constructed values."""

    zero = Angle()
    assert zero.is_valid()
    assert 0.0 == zero.sin().v()
    assert 1.0 == zero.cos().v()
    assert 0.0 == zero.tan()
    assert None == zero.csc()
    assert 1.0 == zero.sec()
    assert None == zero.cot()
    assert 0.0 == zero.to_degrees().v()
    assert 0.0 == zero.to_radians().v()
    assert repr(zero) == 'Angle([ 0.000000, 1.000000 ])'

def test_Angle_y_x_constructor():
    """Test y-x constructed values."""

    one = Angle(1.0, 0.0)
    assert 1.0 == one.sin().v()
    assert 0.0 == one.cos().v()
    assert None == one.tan()
    assert 1.0 == one.csc()
    assert None == one.sec()
    assert 0.0 == one.cot()
    assert Degrees(90.0) == one.to_degrees()
    assert Radians(np.pi / 2) == one.to_radians()
    assert repr(one) == 'Angle([ 1.000000, 0.000000 ])'

    angle_m45 = Angle(-EPSILON, EPSILON)
    assert_almost_equal(-SQRT1_2, angle_m45.sin().v())
    assert_almost_equal(SQRT1_2, angle_m45.cos().v())

def test_Angle_conversion():
    zero = Angle()

    too_small = Angle(-EPSILON / 2.0, EPSILON / 2.0)
    assert too_small.is_valid()
    assert zero == too_small

    small = Angle(Radians(-3.35e7 * EPSILON))
    assert small.is_valid()
    assert -3.35e7 * EPSILON == small.sin().v()
    assert 1.0 == small.cos().v()

    angle_30 = Angle(Radians(PI_3), Radians(PI_6))
    assert angle_30.is_valid()
    assert 0.5 == angle_30.sin().v()
    assert SQRT3 / 2.0 == angle_30.cos().v()
    assert Degrees(30.0) == angle_30.to_degrees()
    assert Radians(PI_6) == angle_30.to_radians()

    angle_45 = Angle(Radians(np.pi / 4.0))
    assert angle_45.is_valid()
    assert SQRT1_2 == angle_45.sin().v()
    assert SQRT1_2 == angle_45.cos().v()
    assert Degrees(45.0) == angle_45.to_degrees()
    assert Radians(np.pi / 4.0) == angle_45.to_radians()

    angle_m45 = Angle(Degrees(-45.0))
    assert angle_m45.is_valid()
    assert -SQRT1_2 == angle_m45.sin().v()
    assert SQRT1_2 == angle_m45.cos().v()
    assert Degrees(-45.0) == angle_m45.to_degrees()
    assert Radians(-np.pi / 4.0) == angle_m45.to_radians()

    angle_60 = Angle(Degrees(-140.0), Degrees(160.0))
    assert angle_30.is_valid()
    assert SQRT3 / 2.0 == angle_60.sin().v()
    assert 0.5 == angle_60.cos().v()
    assert Degrees(60.0) == angle_60.to_degrees()
    # Fails because PI is irrational
    # assert Radians(PI_3) == angle_60.to_radians()
    assert_almost_equal(PI_3, angle_60.to_radians().v())

    angle_d30 = Angle(Degrees(-155.0), Degrees(175.0))
    assert angle_d30.is_valid()
    assert 0.5 == angle_d30.sin().v()
    assert SQRT3 / 2.0 == angle_d30.cos().v()
    assert Degrees(30.0) == angle_d30.to_degrees()
    assert Radians(PI_6) == angle_d30.to_radians()

    angle_120 = Angle(Degrees(120.0))
    assert angle_120.is_valid()
    assert SQRT3 / 2.0 == angle_120.sin().v()
    assert -0.5 == angle_120.cos().v()
    assert Degrees(120.0) == angle_120.to_degrees()
    assert Radians(2.0 * PI_3) == angle_120.to_radians()

    angle_m120 = Angle(Degrees(-120.0))
    assert angle_m120.is_valid()
    assert -SQRT3 / 2.0 == angle_m120.sin().v()
    assert -0.5 == angle_m120.cos().v()
    assert Degrees(-120.0) == angle_m120.to_degrees()
    assert Radians(-2.0 * PI_3) == angle_m120.to_radians()

    angle_m140 = Angle(Degrees(-140.0))
    assert angle_m140.is_valid()
    assert_almost_equal(-0.6427876096865393,  angle_m140.sin().v())
    assert_almost_equal(-0.7660444431189781,  angle_m140.cos().v())
    assert Degrees(-140.0) == angle_m140.to_degrees()

def test_Angle_maths():
    degrees_30 = Angle(Degrees(30.0))
    degrees_60 = Angle(Degrees(60.0))
    degrees_120 = Angle(Degrees(120.0))
    degrees_m120 = -degrees_120

    assert degrees_120 == abs(degrees_m120)
    assert degrees_60 == degrees_m120.opposite()
    assert degrees_120 == degrees_30.quarter_turn_cw()
    assert degrees_30 == degrees_120.quarter_turn_ccw()
    assert degrees_60 == degrees_120.negate_cos()

    assert Degrees(120.0) == (degrees_m120 - degrees_120).to_degrees()
    assert Degrees(-120.0) == (degrees_120 + degrees_120).to_degrees()
    assert Degrees(120.0) == degrees_60.x2().to_degrees()
    assert Degrees(-120.0) == degrees_120.x2().to_degrees()
    assert Degrees(-60.0) == degrees_m120.half().to_degrees()

if __name__ == '__main__':
    pytest.main()
