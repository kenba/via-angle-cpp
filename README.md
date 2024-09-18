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


# Build

Run `cmake` as
```
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 <source dir>
```

to create `compile_commands.json` file for [clangd](https://clangd.llvm.org/).
Then copy `compile_commands.json` back to <source dir> to fix incorrect `clangd` warnings.
