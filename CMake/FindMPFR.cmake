# Try to find the MPFR library
# See http://www.mpfr.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(MPFR 2.3.0)
# to require version 2.3.0 to newer of MPFR.
#
# Once done this will define
#
#  MPFR_FOUND - system has MPFR lib with correct version
#  MPFR_INCLUDE_DIRS - the MPFR include directory
#  MPFR_LIBRARIES - the MPFR library
#  MPFR_VERSION - MPFR version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2010 Jitse Niesen, <jitse@maths.leeds.ac.uk>
# Copyright (c) 2015 Jack Poulson, <jack.poulson@gmail.com>
# Copyright (c) 2020, Mahrud Sayrafi, <mahrud@umn.edu>
# Redistribution and use is allowed according to the terms of the BSD license.

# Set MPFR_FIND_VERSION to 1.0.0 if no minimum version is specified
if(NOT MPFR_FIND_VERSION)
  if(NOT MPFR_FIND_VERSION_MAJOR)
    set(MPFR_FIND_VERSION_MAJOR 1)
  endif()
  if(NOT MPFR_FIND_VERSION_MINOR)
    set(MPFR_FIND_VERSION_MINOR 0)
  endif()
  if(NOT MPFR_FIND_VERSION_PATCH)
    set(MPFR_FIND_VERSION_PATCH 0)
  endif()
  set(MPFR_FIND_VERSION
    "${MPFR_FIND_VERSION_MAJOR}.${MPFR_FIND_VERSION_MINOR}.${MPFR_FIND_VERSION_PATCH}")
endif()

macro(_mpfr_check_version)
  # Query MPFR_VERSION
  file(READ "${MPFR_INCLUDE_DIRS}/mpfr.h" _mpfr_version_header)

  string(REGEX MATCH "define[ \t]+MPFR_VERSION_MAJOR[ \t]+([0-9]+)"
    _mpfr_major_version_match "${_mpfr_version_header}")
  set(MPFR_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+MPFR_VERSION_MINOR[ \t]+([0-9]+)"
    _mpfr_minor_version_match "${_mpfr_version_header}")
  set(MPFR_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+MPFR_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
    _mpfr_patchlevel_version_match "${_mpfr_version_header}")
  set(MPFR_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

  set(MPFR_VERSION
    ${MPFR_MAJOR_VERSION}.${MPFR_MINOR_VERSION}.${MPFR_PATCHLEVEL_VERSION})

  # Check whether found version exceeds minimum required
  if(${MPFR_VERSION} VERSION_LESS ${MPFR_FIND_VERSION})
    set(MPFR_VERSION_OK FALSE)
    message(STATUS "MPFR version ${MPFR_VERSION} found in ${MPFR_INCLUDE_DIRS}, "
                   "but at least version ${MPFR_FIND_VERSION} is required")
  else()
    set(MPFR_VERSION_OK TRUE)
  endif()
endmacro(_mpfr_check_version)

if(NOT MPFR_VERSION_OK)
  set(MPFR_INCLUDE_DIRS NOTFOUND)
  set(MPFR_LIBRARIES NOTFOUND)

  # search first if an MPFRConfig.cmake is available in the system,
  # if successful this would set MPFR_INCLUDE_DIRS and the rest of
  # the script will work as usual
  find_package(MPFR ${MPFR_FIND_VERSION} NO_MODULE QUIET)

  if(NOT MPFR_INCLUDE_DIRS)
    find_path(MPFR_INCLUDE_DIRS NAMES mpfr.h
      HINTS ENV MPFRDIR ENV GMPDIR
      PATHS ${INCLUDE_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/include
      )
  endif()

  if(MPFR_INCLUDE_DIRS)
    _MPFR_check_version()
  endif()

  if(NOT MPFR_LIBRARIES)
    find_library(MPFR_LIBRARIES NAMES mpfr
      HINTS ENV MPFRDIR ENV GMPDIR
      PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
      )
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_INCLUDE_DIRS MPFR_LIBRARIES MPFR_VERSION_OK)

  mark_as_advanced(MPFR_INCLUDE_DIRS MPFR_LIBRARIES)

endif()
