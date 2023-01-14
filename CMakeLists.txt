cmake_minimum_required(VERSION 3.10)

project(cgshop2023 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED true)

set(SRC ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
include_directories("${SRC}")

#set(CMAKE_CXX_FLAGS "-Wall -pthread -O3 -mavx -fopenmp -g -fsanitize=address")
set(CMAKE_CXX_FLAGS "-Wall -pthread -O3 -mavx -fopenmp -g")

include(./cmake/CPM.cmake)  # Package manager for simple requirements.

# run `conan install .. --build=missing` in `./cmake`
list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

CPMAddPackage("gh:fmtlib/fmt#9.1.0")  # fmt for nice strings

find_package(Threads)
#find_package(X11)
find_package(CGAL REQUIRED COMPONENTS Core)
#find_package(CGAL REQUIRED COMPONENTS Core Qt5)
#find_package(Qt5 QUIET COMPONENTS Qt5 Core)

add_definitions("-DCGAL_USE_BASIC_VIEWER")

file(GLOB SOURCES_simple ${SRC}/*.cpp ${SRC}/cgshop2023_core/*.cpp)

add_executable(simple ${SOURCES_simple})

#target_link_libraries(simple PUBLIC CGAL::CGAL CGAL::CGAL_Qt5 fmt::fmt)
target_link_libraries(simple PUBLIC CGAL::CGAL fmt::fmt)
