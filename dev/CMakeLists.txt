cmake_minimum_required(VERSION 3.10)
project(dev)

set(CMAKE_CXX_STANDARD 14)

add_library(dev SHARED library.cpp library.h Grid.cpp Grid.h)