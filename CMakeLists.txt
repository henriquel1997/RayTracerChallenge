cmake_minimum_required(VERSION 3.13)
project(RayTracerChallenge)

set(CMAKE_CXX_STANDARD 11)

add_executable(RayTracerChallenge main.cpp tuples.h util.h projectile.h color.h canvas.h matrix.h transformations.h clock.h)