cmake_minimum_required(VERSION 3.12)
project(PracticeCode)

set(CMAKE_CXX_STANDARD 14)

add_executable(PracticeCode main.cpp Grid.cpp Grid.h Wavefunction.cpp Wavefunction.h Potential.cpp Potential.h)

find_package(GSL REQUIRED)
target_link_libraries(PracticeCode GSL::gsl GSL::gslcblas)