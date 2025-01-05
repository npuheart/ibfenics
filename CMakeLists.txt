cmake_minimum_required(VERSION 3.4...3.18)
project(cpp)

# Find pybind11
find_package(pybind11 REQUIRED)

# Find dolfin
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

add_definitions("-Wall -g")
add_definitions("-O3")

# pybind11_add_module(cpp src/main.cpp)
pybind11_add_module(cpp src/IB.cpp src/IBMesh.cpp src/IBInterpolation.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(cpp
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})


add_subdirectory(src/gauss_quadrature)

# Target the library
target_link_libraries(cpp PRIVATE dolfin gauss_quadrature)
