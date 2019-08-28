cmake_minimum_required(VERSION 3.10)

project(QuantumToolbox)

set(CMAKE_CXX_STANDARD 14)
#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
#if (COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#else()
#    CHECK_CXX_COMPILER_FLAG("/arch:AVX" COMPILER_OPT_ARCH_AVX_SUPPORTED)
#    if(COMPILER_OPT_ARCH_AVX_SUPPORTED)
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
#    endif()
#endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O5 -msse -msse2 -msse3")
set(CMAKE_CXX_FLAGS " -g -O0" )

add_executable(QuantumToolbox main.cpp Grid.cpp Grid.h Wavefunction.cpp Wavefunction.h Potential.cpp Potential.h Extras.cpp Extras.h System.cpp System.h Plotter.cpp Plotter.h)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS system regex iostreams)

#include_directories ($ENV{EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

#target_link_libraries(QuantumToolbox Eigen3::Eigen GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(QuantumToolbox GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})