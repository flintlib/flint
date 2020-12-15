# Try to find BLAS, including cblas.h
#
# Once done this will define
#  CBLAS_FOUND        - system has a BLAS library
#  CBLAS_INCLUDE_DIRS - the header directory containing cblas.h
#  CBLAS_LIBRARIES    - the CBLAS library

# Copyright (c) 2020, Mahrud Sayrafi, <mahrud@umn.edu>
# Redistribution and use is allowed according to the terms of the BSD license.

find_path(CBLAS_INCLUDE_DIRS NAMES cblas.h
  HINTS CBLAS_ROOT ENV CBLAS_ROOT
  PATHS ${INCLUDE_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/include
  PATH_SUFFIXES openblas cblas blis
  )

find_library(CBLAS_LIBRARIES NAMES accelerate openblas cblas blas blis
  HINTS CBLAS_ROOT ENV CBLAS_ROOT
  PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
  PATH_SUFFIXES openblas cblas blis
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS
  "Could NOT find a BLAS compatible library or 'cblas.h', install BLAS or set CBLAS_ROOT."
  CBLAS_INCLUDE_DIRS CBLAS_LIBRARIES)

mark_as_advanced(CBLAS_LIBRARIES CBLAS_INCLUDE_DIRS)
