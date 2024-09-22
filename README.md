# via-angle-cpp

A C++ library for performing accurate and efficient trigonometry calculations.

## Description

The standard trigonometry functions: `sin`, `cos`, `tan`, etc.
[give unexpected results for well-known angles](https://stackoverflow.com/questions/31502120/sin-and-cos-give-unexpected-results-for-well-known-angles#answer-31525208).
This is because the functions use parameters with `radians` units instead of `degrees`.
The conversion from `degrees` to `radians` suffers from
[round-off error](https://en.wikipedia.org/wiki/Round-off_error) due to
`radians` being based on the irrational number π.

This library provides a [sincos](include/via/trig.hpp#sincos) function to calculate more
accurate values than the standard `sin` and `cos` functions for angles in radians
and a [sincosd](include/via/trig.hpp#sincosd) function to calculate more accurate values
for angles in degrees.

The library also provides an [Angle](#angle) class which represents an angle
by its sine and cosine as the coordinates of a
[unit circle](https://en.wikipedia.org/wiki/Unit_circle), see *Figure 1*.

![Unit circle](https://upload.wikimedia.org/wikipedia/commons/thumb/7/72/Sinus_und_Kosinus_am_Einheitskreis_1.svg/250px-Sinus_und_Kosinus_am_Einheitskreis_1.svg.png)
*Figure 1 Unit circle formed by cos *θ* and sin *θ**

The `Angle` class enables more accurate calculations of angle rotations and
conversions to and from `degrees` or `radians`.

## Features

* `Degrees`, `Radians` and `Angle` types;
* functions for accurately calculating sines and cosines of angles in `Degrees` or `Radians`
using [remquo](https://pubs.opengroup.org/onlinepubs/9699919799/functions/remquo.html);
* functions for accurately calculating sines and cosines of differences of angles in `Degrees` or `Radians`
using the [2Sum](https://en.wikipedia.org/wiki/2Sum) algorithm;
* functions for accurately calculating sums and differences of `Angles` using
[trigonometric identities](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Angle_sum_and_difference_identities);
* and some [spherical trigonometry](https://en.wikipedia.org/wiki/Spherical_trigonometry) functions.

The following example shows the `round-off error` inherent in calculating angles in `radians`.  
It calculates the correct sine and cosine for 60° and converts them back
precisely to 60°, but it fails to convert them to the precise angle in `radians`: π/3.

```C++
#include "via/angle.hpp"
#include <boost/test/unit_test.hpp>

using namespace via;

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
constexpr auto CALCULATION_TOLERANCE{101 * EPSILON};
} // namespace

BOOST_AUTO_TEST_SUITE(Test_angle)

BOOST_AUTO_TEST_CASE(test_Angle_conversion) {
  const Angle angle_60(Degrees(60.0));
  BOOST_CHECK(angle_60.is_valid());
  BOOST_CHECK_EQUAL(trig::SQRT3<double> / 2, angle_60.sin().v());
  BOOST_CHECK_EQUAL(0.5, angle_60.cos().v());
  BOOST_CHECK_EQUAL(Degrees(60.0), angle_60.to_degrees());
  // Fails because PI is irrational
  // BOOST_CHECK_EQUAL(Radians(trig::PI_3<double>), angle_60.to_radians());
  BOOST_CHECK_CLOSE(trig::PI_3<double>, angle_60.to_radians().v(),
                    CALCULATION_TOLERANCE);
}

BOOST_AUTO_TEST_SUITE_END()
```
The following example calculates the sine and cosine between the difference
of two angles in `degrees`: -155° - 175°.  
It is more accurate than calling the `Angle` `From` trait in the example above
with the difference in `degrees`.  
It is particularly useful for implementing the
[Haversine formula](https://en.wikipedia.org/wiki/Haversine_formula)
which requires sines and cosines of both longitude and latitude differences.  
Note: in this example sine and cosine of 30° are converted precisely to π/6.

```C++
#include "via/angle.hpp"
#include <boost/test/unit_test.hpp>

using namespace via;

namespace {
constexpr auto EPSILON{std::numeric_limits<double>::epsilon()};
constexpr auto CALCULATION_TOLERANCE{101 * EPSILON};
} // namespace

BOOST_AUTO_TEST_SUITE(Test_angle_difference)

BOOST_AUTO_TEST_CASE(test_Angle_difference_conversion) {
  const Angle angle_d30(Degrees(-155.0), Degrees(175.0));
  BOOST_CHECK(angle_d30.is_valid());
  BOOST_CHECK_EQUAL(0.5, angle_d30.sin().v());
  BOOST_CHECK_EQUAL(trig::SQRT3<double> / 2, angle_d30.cos().v());
  BOOST_CHECK_EQUAL(Degrees(30.0), angle_d30.to_degrees());
  BOOST_CHECK_EQUAL(Radians(trig::PI_6<double>), angle_d30.to_radians());
}

BOOST_AUTO_TEST_SUITE_END()
```

## Design

### Trigonometry Functions

The [trig](include/via/trig.hpp) namespace contains accurate and efficient trigonometry functions.

### Angle

The `Angle` struct represents an angle by its sine and cosine instead of in
`degrees` or `radians`, see *Figure 2*.  

![Angle Class Diagram](docs/images/angle_class_diagram.svg)  
*Figure 2 Angle Class Diagram*

This representation an angle makes functions such as
rotating an angle +/-90° around the unit circle or calculating the opposite angle;
simple, accurate and efficient since they just involve changing the signs
and/or positions of the `sin` and `cos` values.

`Angle` `Add` and `Sub` traits are implemented using
[angle sum and difference](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Angle_sum_and_difference_identities)
trigonometric identities, 
while `Angle` [double](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Double-angle_formulae)
and [half](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Half-angle_formulae) methods use other
trigonometric identities.

The `sin` and `cos` fields of `Angle` are `UnitNegRange`s:,
a [newtype](https://rust-unofficial.github.io/patterns/patterns/behavioural/newtype.html)
with values in the range -1.0 to +1.0 inclusive.  

# Build

Run `cmake` as
```
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 <source dir>
```

to create `compile_commands.json` file for [clangd](https://clangd.llvm.org/).  
Then copy `compile_commands.json` back to <source dir> to fix incorrect `clangd` warnings.

## License

`angle-rs` is provided under a BOOST license, see [LICENSE](LICENSE).
