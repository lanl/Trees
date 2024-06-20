cmake_minimum_required(VERSION 3.5)

# Project Information
project(DUET)
set(DUET_VERSION_MAJOR 1)
set(DUET_VERSION_MINOR 1.0)
set(DUET_VERSION "${DUET_VERSION_MAJOR}.${DUET_VERSION_MINOR}")

# Build
set(DUET_SRCS MakeFFInputs_DUET.f90 variables_DUET.f90
  define_variables_DUET.f90 io_DUET.f90 full_DUET.f90)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../opt)
add_library(DUET STATIC ${DUET_SRCS})
