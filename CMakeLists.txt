cmake_minimum_required(VERSION 2.8.12)

project(QuantumToolbox)

set(CMAKE_CXX_STANDARD 17)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_COMPILER g++)

CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
if (COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O5 -msse -msse2 -msse3 -mmmx -m3dnow")
else()
    CHECK_CXX_COMPILER_FLAG("/arch:AVX" COMPILER_OPT_ARCH_AVX_SUPPORTED)
    if(COMPILER_OPT_ARCH_AVX_SUPPORTED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX -O5 -msse -msse2 -msse3 -mmmx -m3dnow")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O5 -msse -msse2 -msse3 -mmmx -m3dnow")
    endif()
endif()

#set(CMAKE_CXX_FLAGS " -g -O0" )

add_executable(${PROJECT_NAME} main.cpp Grid.cpp Grid.h Wavefunction.cpp Wavefunction.h Potential.cpp Potential.h Extras.cpp Extras.h System.cpp System.h Plotter.cpp Plotter.h)
#add_library(${PROJECT_NAME} Grid.cpp Grid.h Wavefunction.cpp Wavefunction.h Potential.cpp Potential.h Extras.cpp Extras.h System.cpp System.h Plotter.cpp Plotter.h)

#include_directories ($ENV{EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories($ENV{GSL_INCLUDE_DIRS})
include_directories(/eigen)

#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS system regex iostreams filesystem)
find_package(OpenMP REQUIRED) # Find the package

#target_link_libraries(${PROJECT_NAME} Eigen3::Eigen GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX GSL::gsl GSL::gslcblas ${Boost_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX ${Boost_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

