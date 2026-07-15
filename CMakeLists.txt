cmake_minimum_required(VERSION 3.10)
project(rm_synthesis Fortran)

# Set Fortran standard and flags
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -fbounds-check")

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags for gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=gnu -fallow-argument-mismatch -ffree-line-length-none")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -assume byterecl")
endif()

# Find CFITSIO library
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(CFITSIO QUIET cfitsio)
endif()

if(NOT CFITSIO_FOUND)
  # Manual search
  find_library(CFITSIO_LIBRARY cfitsio)
  find_path(CFITSIO_INCLUDE_DIR fitsio.h)
  
  if(CFITSIO_LIBRARY AND CFITSIO_INCLUDE_DIR)
    set(CFITSIO_FOUND TRUE)
    set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY})
    set(CFITSIO_INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})
  endif()
endif()

if(NOT CFITSIO_FOUND)
  message(FATAL_ERROR "CFITSIO library not found. Install it with: sudo apt-get install libcfitsio-dev")
endif()

message(STATUS "CFITSIO found: ${CFITSIO_LIBRARIES}")

# Source files
set(SOURCES
  src/rm_synthesis_mod.f90
  src/rm_synthesis.f90
  src/myfits_info.f
  src/printerror.f
)

# Create executable
add_executable(rm_synthesis ${SOURCES})

# Link CFITSIO
target_link_libraries(rm_synthesis PRIVATE ${CFITSIO_LIBRARIES})
if(CFITSIO_INCLUDE_DIRS)
  target_include_directories(rm_synthesis PRIVATE ${CFITSIO_INCLUDE_DIRS})
endif()

# Installation
install(TARGETS rm_synthesis DESTINATION bin)
install(DIRECTORY cfg/ DESTINATION share/rm_synthesis/cfg)

# Print build info
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
