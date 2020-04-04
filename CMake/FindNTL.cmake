# Try to find the NTL libraries
# See https://www.shoup.net/ntl/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(NTL 10.5.0)
# to require version 10.5.0 to newer of NTL.
#
# Once done this will define
#
#  NTL_FOUND		- system has the NTL library with correct version
#  NTL_INCLUDE_DIR	- the NTL include directory
#  NTL_LIBRARIES	- the NTL library
#  NTL_VERSION		- NTL version

# Copyright (c) 2020, Mahrud Sayrafi, <mahrud@umn.edu>
#
# Redistribution and use is allowed according to the terms of the BSD license.

# Set NTL_FIND_VERSION to 1.0.0 if no minimum version is specified
if(NOT NTL_FIND_VERSION)
  if(NOT NTL_FIND_VERSION_MAJOR)
    set(NTL_FIND_VERSION_MAJOR 1)
  endif()
  if(NOT NTL_FIND_VERSION_MINOR)
    set(NTL_FIND_VERSION_MINOR 0)
  endif()
  if(NOT NTL_FIND_VERSION_PATCH)
    set(NTL_FIND_VERSION_PATCH 0)
  endif()

  set(NTL_FIND_VERSION "${NTL_FIND_VERSION_MAJOR}.${NTL_FIND_VERSION_MINOR}.${NTL_FIND_VERSION_PATCH}")
endif()

macro(_NTL_check_version)
  # Query NTL_VERSION
  file(READ "${NTL_INCLUDE_DIR}/NTL/version.h" _NTL_version_header)

  string(REGEX MATCH "define[ \t]+NTL_MAJOR_VERSION[ \t]+\\(([0-9]+)\\)"
    _NTL_major_version_match "${_NTL_version_header}")
  set(NTL_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+NTL_MINOR_VERSION[ \t]+\\(([0-9]+)\\)"
    _NTL_minor_version_match "${_NTL_version_header}")
  set(NTL_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+NTL_REVISION[ \t]+\\(([0-9]+)\\)"
    _NTL_patchlevel_version_match "${_NTL_version_header}")
  set(NTL_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

  set(NTL_VERSION
    ${NTL_MAJOR_VERSION}.${NTL_MINOR_VERSION}.${NTL_PATCHLEVEL_VERSION})

  # Check whether found version exceeds minimum required
  if(${NTL_VERSION} VERSION_LESS ${NTL_FIND_VERSION})
    set(NTL_VERSION_OK FALSE)
    message(STATUS "NTL version ${NTL_VERSION} found in ${NTL_INCLUDE_DIR}, "
                   "but at least version ${NTL_FIND_VERSION} is required")
  else()
    set(NTL_VERSION_OK TRUE)
  endif()
endmacro(_NTL_check_version)

if(NOT NTL_FOUND)
  set(NTL_INCLUDE_DIR NOTFOUND)
  set(NTL_LIBRARIES NOTFOUND)

  # search first if an NTLConfig.cmake is available in the system,
  # if successful this would set NTL_INCLUDE_DIR and the rest of
  # the script will work as usual
  find_package(NTL ${NTL_FIND_VERSION} NO_MODULE QUIET)

  if(NOT NTL_INCLUDE_DIR)
    find_path(NTL_INCLUDE_DIR NAMES NTL/version.h
      HINTS ENV NTLDIR
      PATHS ${INCLUDE_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/include
      PATH_SUFFIXES NTL
      )
  endif()

  if(NTL_INCLUDE_DIR)
    _NTL_check_version()
  endif()

  if(NOT NTL_LIBRARIES)
    find_library(NTL_LIBRARIES NAMES ntl
      HINTS ENV NTLDIR
      PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
      )
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(NTL DEFAULT_MSG NTL_INCLUDE_DIR NTL_LIBRARIES NTL_VERSION_OK)

  mark_as_advanced(NTL_INCLUDE_DIR NTL_LIBRARIES)

endif()
