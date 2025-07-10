cmake_minimum_required(VERSION 3.5)

# Project Information
project(TREES)
set(TREES_VERSION_MAJOR 1)

# Compiler information
enable_language(Fortran)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set(FF_FLAGS "-O2 -g -ffixed-line-length-none")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set(FF_FLAGS "-O2 -g -extend-source")
elseif (Fortran_COMPILER_NAME MATCHES "ftn.*")
  set(FF_FLAGS "-O2 -g -extend_source")
endif()

#find_package(MPI REQUIRED)
#set(-Wall -Wextra -Werror -fcheck=all) # Flags for debugging

# Module dependecies
add_subdirectory(ModularDUET)

# Build
set(Trees_SRCS
  variables.F90 main.F90 define_variables.F90
  power_parser.F90 io.F90 metryc.F90 shapes.F90
  distribution.F90 fuels_create.F90
  fuel_read.F90 trees_read.F90)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/opt)
add_executable(trees.exe ${Trees_SRCS})
set_target_properties(trees.exe PROPERTIES COMPILE_FLAGS "${FF_FLAGS}")
# target_compile_options(firerad PRIVATE -Wall -Wextra)
target_link_libraries(trees.exe DUET)

# Make Clean Instructions
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
  "opt/*;trees.exe;DUET/libDUET.a")
