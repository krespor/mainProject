cmake_minimum_required(VERSION 3.13)
project(NewFEM)

set(CMAKE_CXX_STANDARD 14)

add_executable(NewFEM main.cpp structurs/arguments.h structurs/mesh.h structurs/paths.h Solver.cpp Solver.h structurs/localMatrix.h structurs/methods.h SUPG.cpp SUPG.h NavierStokes.cpp NavierStokes.h)