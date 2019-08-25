cmake_minimum_required(VERSION 3.10)

project(QuantumToolbox)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g -O0" )

#add_executable(QuantumToolbox main.cpp Grid.cpp Grid.h Wavefunction.cpp Wavefunction.h Potential.cpp Potential.h Extras.cpp Extras.h System.cpp System.h Plotter.cpp Plotter.h)
add_executable(QuantumToolbox main.cpp Grid.cpp Grid.h Wavefunction.cpp Wavefunction.h Potential.cpp Potential.h Extras.cpp Extras.h )
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS system regex iostreams)

include_directories ($ENV{EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

#target_link_libraries(QuantumToolbox Eigen3::Eigen GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(QuantumToolbox GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})