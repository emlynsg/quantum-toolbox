cmake_minimum_required(VERSION 3.12)
project(QuantumToolbox)

set(CMAKE_CXX_STANDARD 14)

add_executable(QuantumToolbox main.cpp Grid.cpp Grid.h Wavefunction.cpp Wavefunction.h Potential.cpp Potential.h Extras.cpp Extras.h System.cpp System.h Plotter.cpp Plotter.h)

find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS system regex iostreams)
include_directories(${Boost_INCLUDE_DIRS})

target_link_libraries(QuantumToolbox GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})