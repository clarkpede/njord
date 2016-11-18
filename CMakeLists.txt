##
#  CMake script for the backward facing step program
##

cmake_minimum_required(VERSION 2.6)

project(NJORD C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH}})

# PETSc dependencies
find_package(PETSc REQUIRED)
message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")

FILE(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/*.c ${CMAKE_SOURCE_DIR}/src/*.cc)
FILE(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/*.h)

# Essential: include our directories first otherwise we can get internal headers from some installed path
include_directories (${CMAKE_SOURCE_DIR}/include ${PETSC_INCLUDES})
add_definitions (${PETSC_DEFINITIONS})

add_executable(njord ${SOURCES} ${HEADERS})

target_link_libraries(njord ${PETSC_LIBRARIES})



