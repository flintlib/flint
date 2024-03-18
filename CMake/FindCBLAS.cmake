# Try to find BLAS, including cblas.h
#
# Once done this will define
#  CBLAS_FOUND        - system has a BLAS library
#  CBLAS_INCLUDE_DIR  - the header directory containing cblas.h
#  CBLAS_LIBRARY      - the CBLAS library

# Copyright (c) 2020, Mahrud Sayrafi, <mahrud@umn.edu>
# Redistribution and use is allowed according to the terms of the BSD license.

find_path(CBLAS_INCLUDE_DIR NAMES cblas.h
  HINTS CBLAS_ROOT ENV CBLAS_ROOT
  PATHS ${INCLUDE_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/include
  PATH_SUFFIXES openblas cblas blis flexiblas
  )

find_library(CBLAS_LIBRARY NAMES accelerate openblas cblas blas blis flexiblas
  HINTS CBLAS_ROOT ENV CBLAS_ROOT
  PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
  PATH_SUFFIXES openblas cblas blis flexiblas
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( CBLAS
  FOUND_VAR CBLAS_FOUND
  REQUIRED_VARS
    CBLAS_LIBRARY
    CBLAS_INCLUDE_DIR
  )

if(CBLAS_FOUND)
  set(CBLAS_INCLUDE_DIRS ${CBLAS_INCLUDE_DIR})
  set(CBLAS_LIBRARIES ${CBLAS_LIBRARY})
  if(NOT TARGET CBLAS::CBLAS)
    add_library(CBLAS::CBLAS UNKNOWN IMPORTED)
    set_target_properties( CBLAS::CBLAS
      PROPERTIES
        IMPORTED_LOCATION "${CBLAS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CBLAS_INCLUDE_DIR}" 
      )
  endif()
  mark_as_advanced(CBLAS_ROOT)
endif()

mark_as_advanced(CBLAS_LIBRARY CBLAS_INCLUDE_DIR)
