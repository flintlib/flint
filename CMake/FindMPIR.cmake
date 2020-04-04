# Try to find the MPIR libraries
# See http://www.mpir.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(MPIR 3.0.0)
# to require version 3.0.0 to newer of MPIR.
#
# Once done this will define
#
#  MPIR_FOUND             - system has MPIR lib
#  MPIR_INCLUDE_DIRS      - the MPIR include directory
#  MPIR_LIBRARIES         - Libraries needed to use MPIR
#  MPIR_VERSION           - MPIR version

# Copyright (c) 2020, Mahrud Sayrafi, <mahrud@umn.edu>
# Redistribution and use is allowed according to the terms of the BSD license.

# Set MPIR_FIND_VERSION to 1.0.0 if no minimum version is specified
if(NOT MPIR_FIND_VERSION)
  if(NOT MPIR_FIND_VERSION_MAJOR)
    set(MPIR_FIND_VERSION_MAJOR 1)
  endif()
  if(NOT MPIR_FIND_VERSION_MINOR)
    set(MPIR_FIND_VERSION_MINOR 0)
  endif()
  if(NOT MPIR_FIND_VERSION_PATCH)
    set(MPIR_FIND_VERSION_PATCH 0)
  endif()
  set(MPIR_FIND_VERSION
    "${MPIR_FIND_VERSION_MAJOR}.${MPIR_FIND_VERSION_MINOR}.${MPIR_FIND_VERSION_PATCH}")
endif()

macro(_mpir_check_version)
  # Query MPIR_VERSION
  file(READ "${MPIR_INCLUDE_DIRS}/mpir.h" _mpir_version_header)

  string(REGEX MATCH "define[ \t]+__MPIR_VERSION[ \t]+([0-9]+)"
    _mpir_major_version_match "${_mpir_version_header}")
  set(MPIR_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+__MPIR_VERSION_MINOR[ \t]+([0-9]+)"
    _mpir_minor_version_match "${_mpir_version_header}")
  set(MPIR_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+__MPIR_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
    _mpir_patchlevel_version_match "${_mpir_version_header}")
  set(MPIR_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

  set(MPIR_VERSION
    ${MPIR_MAJOR_VERSION}.${MPIR_MINOR_VERSION}.${MPIR_PATCHLEVEL_VERSION})

  # Check whether found version exceeds minimum required
  if(${MPIR_VERSION} VERSION_LESS ${MPIR_FIND_VERSION})
    set(MPIR_VERSION_OK FALSE)
    message(STATUS "MPIR version ${MPIR_VERSION} found in ${MPIR_INCLUDE_DIRS}, "
                   "but at least version ${MPIR_FIND_VERSION} is required")
  else()
    set(MPIR_VERSION_OK TRUE)
  endif()
endmacro(_mpir_check_version)

if(NOT MPIR_VERSION_OK)
  set(MPIR_INCLUDE_DIRS NOTFOUND)
  set(MPIR_LIBRARIES NOTFOUND)

  # search first if an MPIRConfig.cmake is available in the system,
  # if successful this would set MPIR_INCLUDE_DIRS and the rest of
  # the script will work as usual
  find_package(MPIR ${MPIR_FIND_VERSION} NO_MODULE QUIET)

  if(NOT MPIR_INCLUDE_DIRS)
    find_path(MPIR_INCLUDE_DIRS NAMES mpir.h
      HINTS ENV MPIRDIR ENV GMPDIR
      PATHS ${INCLUDE_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/include
      )
  endif()

  if(MPIR_INCLUDE_DIRS)
    _MPIR_check_version()
  endif()

  if(NOT MPIR_LIBRARIES)
    find_library(MPIR_LIBRARIES NAMES mpir
      HINTS ENV MPIRDIR ENV GMPDIR
      PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
      )
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MPIR DEFAULT_MSG MPIR_INCLUDE_DIRS MPIR_LIBRARIES MPIR_VERSION_OK)

  mark_as_advanced(MPIR_INCLUDE_DIRS MPIR_LIBRARIES)

endif()
