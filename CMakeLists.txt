cmake_minimum_required(VERSION 3.15)
project(2d_calc C)

set(CMAKE_C_STANDARD 99)

add_executable(2d_calc main.c)

target_link_libraries(2d_calc m)