cmake_minimum_required (VERSION 3.22.1)
project (MeshGen)

add_subdirectory(extern/pybind11)
pybind11_add_module(${PROJECT_NAME} SHARED src/subdivision/subdivide.cpp src/subdivision/subdivide.h)
