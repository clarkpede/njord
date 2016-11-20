##
#  CMake script for the backward facing step program
##

cmake_minimum_required(VERSION 2.6)

project(NJORD C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH}})

# PETSc dependencies
find_package(PETSc REQUIRED)
find_package(MPI REQUIRED)
message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")

# BOOST Test
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories($(Boost_INCLUDE_DIRS), ${MPI_INCLUDE_PATH})

# Source files
FILE(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/*.c ${CMAKE_SOURCE_DIR}/src/*.cc)
FILE(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/*.h)

# Test files
enable_testing()
FILE(GLOB TEST_SOURCES ${CMAKE_SOURCE_DIR}/tests/*.cc)

# Run through each source
add_executable(unit_tests ${TEST_SOURCES})
target_link_libraries(unit_tests ${Boost_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_C_LIBRARIES})
# Loop all the test suites and run each individually
foreach(testSrc ${TEST_SOURCES})
    IF (NOT ${testSrc} MATCHES "main_test.cc") 
        get_filename_component(testName ${testSrc} NAME_WE)
        add_test(NAME ${testName} COMMAND ./unit_tests -run_test=${testName})
    ENDIF (NOT ${testSrc} MATCHES "main_test.cc") 
endforeach(testSrc)


# Essential: include our directories first otherwise we can get internal headers from some installed path
include_directories (${CMAKE_SOURCE_DIR}/include ${PETSC_INCLUDES})
add_definitions (${PETSC_DEFINITIONS})

# Build and link
add_executable(njord ${SOURCES} ${HEADERS})
target_link_libraries(njord ${PETSC_LIBRARIES})



