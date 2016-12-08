##
#  CMake script for the backward facing step program
##

cmake_minimum_required(VERSION 3.1)

project(NJORD C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Intel compiler support is only available on CMake 3.5 or later
# This is a hack to get around the lack of compatibility
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH}})

# PETSc dependencies
find_package(PETSc REQUIRED)
find_package(MPI REQUIRED)
message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")

# BOOST Test
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories($(Boost_INCLUDE_DIRS), ${MPI_INCLUDE_PATH})

# Source files
set(MAIN "${CMAKE_SOURCE_DIR}/src/main.cc")
file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/*.c ${CMAKE_SOURCE_DIR}/src/*.cc)
file(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/*.h)
# Separate main from the source files listing
list(REMOVE_ITEM SOURCES ${CMAKE_SOURCE_DIR}/src/main.cc)

# Test files
enable_testing()
file(GLOB TEST_SOURCES ${CMAKE_SOURCE_DIR}/tests/*.cc)

# Essential: include our directories first otherwise we can get internal headers from some installed path
include_directories (${CMAKE_SOURCE_DIR}/include ${PETSC_INCLUDES})
add_definitions (${PETSC_DEFINITIONS})

# Run through each source
add_executable(unit_tests ${TEST_SOURCES} ${SOURCES} ${HEADERS})
target_link_libraries(unit_tests ${Boost_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_C_LIBRARIES})
# Loop all the test suites and run each individually
foreach(testSrc ${TEST_SOURCES})
    if (NOT ${testSrc} MATCHES "main_test.cc") 
        get_filename_component(testName ${testSrc} NAME_WE)
        add_test(NAME ${testName} COMMAND ./unit_tests -run_test=${testName})
    endif (NOT ${testSrc} MATCHES "main_test.cc") 
endforeach(testSrc)

# Build and link
add_executable(njord ${MAIN} ${SOURCES} ${HEADERS})
target_link_libraries(njord ${PETSC_LIBRARIES})



