cmake_minimum_required(VERSION 3.15)
project(2d_calc C)

set(CMAKE_C_STANDARD 99)

add_executable(2d_calc src/main.c src/include/global.h src/test.c src/basis.c src/include/basis.h src/utils.c src/include/chemistry.h src/chemistry.c src/include/physics.h src/physics.c src/include/parameters.h src/parameters.c)

target_link_libraries(2d_calc m)