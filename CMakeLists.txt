# Copyright (c) 2024 Ken Barker. All Rights Reserved.
# (ken dot barker at via-technology dot co dot uk)
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required (VERSION 3.30)
cmake_policy(SET CMP0074 NEW)
project (via-angle-cpp)

option(CODE_COVERAGE "Add gcc code coverage options." OFF)

add_library(${PROJECT_NAME} INTERFACE)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_23)

# define _USE_MATH_DEFINES to define M_PI etc.
target_compile_definitions(${PROJECT_NAME} INTERFACE _USE_MATH_DEFINES)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

# Add the C++ Guidelines Support Library include files
set(CPP_GSL_DIR $ENV{CPP_GSL_DIR})
include_directories(${CPP_GSL_DIR}/include)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)
if(Boost_FOUND)
  add_executable(${PROJECT_NAME}_test "")

  file(GLOB_RECURSE INCLUDE_FILES include/via/*.hpp)
  target_sources(${PROJECT_NAME}_test
    PRIVATE
      ${INCLUDE_FILES}
  )

  target_sources(${PROJECT_NAME}_test
    PRIVATE
      tests/test_main.cpp

      tests/test_two_sum.cpp
      tests/test_trig_double.cpp
      tests/test_angle_double.cpp
  )

  target_compile_definitions(${PROJECT_NAME}_test PRIVATE BOOST_TEST_DYN_LINK)
  target_include_directories(${PROJECT_NAME}_test PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME}_test
                        PRIVATE
                        ${PROJECT_NAME}
                        Boost::unit_test_framework)

  if (MSVC)
    target_compile_options(${PROJECT_NAME}_test PRIVATE /W4)
  else()
    target_compile_options(${PROJECT_NAME}_test PRIVATE -Wall -Wextra -Wpedantic)

    if (CODE_COVERAGE)
      target_compile_options(${PROJECT_NAME}_test PRIVATE --coverage)
      target_link_libraries(${PROJECT_NAME}_test PRIVATE --coverage)
    endif()
  endif()

  add_test(NAME via_angle.test COMMAND  ${PROJECT_NAME}_test)
endif(Boost_FOUND)
