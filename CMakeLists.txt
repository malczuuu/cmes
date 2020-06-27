cmake_minimum_required(VERSION 3.10)
project(cmes)

set(CMAKE_CXX_STANDARD 14)

add_executable(cmes
        src/func1d.hpp
        src/func1d.cpp
        src/solver.hpp
        src/solver.cpp
        src/main.hpp
        src/main.cpp
        )
