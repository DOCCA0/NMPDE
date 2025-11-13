cmake_minimum_required(VERSION 3.12.0)
project(lab01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(stokes-time src/stokes-time.cpp src/main-stokes-time.cpp)
deal_ii_setup_target(stokes-time)
