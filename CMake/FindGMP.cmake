# Try to find the GNU Multiple Precision Arithmetic Library (GMP)
# See http://gmplib.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(GMP 6.1.0)
# to require version 6.1.0 to newer of GMP.
#
# Once done this will define
#
#  GMP_FOUND             - system has GMP lib
#  GMP_INCLUDE_DIRS      - the GMP include directory
#  GMP_LIBRARIES         - Libraries needed to use GMP
#  GMP_VERSION           - GMP version

# Copyright (c) 2016 Jack Poulson, <jack.poulson@gmail.com>
# Copyright (c) 2020, Mahrud Sayrafi, <mahrud@umn.edu>
# Redistribution and use is allowed according to the terms of the BSD license.

# Set GMP_FIND_VERSION to 1.0.0 if no minimum version is specified
if(NOT GMP_FIND_VERSION)
  if(NOT GMP_FIND_VERSION_MAJOR)
    set(GMP_FIND_VERSION_MAJOR 5)
  endif()
  if(NOT GMP_FIND_VERSION_MINOR)
    set(GMP_FIND_VERSION_MINOR 1)
  endif()
  if(NOT GMP_FIND_VERSION_PATCH)
    set(GMP_FIND_VERSION_PATCH 0)
  endif()
  set(GMP_FIND_VERSION
    "${GMP_FIND_VERSION_MAJOR}.${GMP_FIND_VERSION_MINOR}.${GMP_FIND_VERSION_PATCH}")
endif()

macro(_gmp_check_version)
  # Since the GMP version macros may be in a file included by gmp.h of the form
  # gmp-.*[_]?.*.h (e.g., gmp-x86_64.h), we search each of them.
  file(GLOB GMP_HEADERS "${GMP_INCLUDE_DIRS}/gmp.h" "${GMP_INCLUDE_DIRS}/gmp-*.h")
  foreach(gmp_header_filename ${GMP_HEADERS})
    file(READ "${gmp_header_filename}" _gmp_version_header)
    string(REGEX MATCH
      "define[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+)" _gmp_major_version_match
      "${_gmp_version_header}")
    if(_gmp_major_version_match)
      set(GMP_MAJOR_VERSION "${CMAKE_MATCH_1}")
      string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)"
        _gmp_minor_version_match "${_gmp_version_header}")
      set(GMP_MINOR_VERSION "${CMAKE_MATCH_1}")
      string(REGEX MATCH "define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
        _gmp_patchlevel_version_match "${_gmp_version_header}")
      set(GMP_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")
      set(GMP_VERSION
        ${GMP_MAJOR_VERSION}.${GMP_MINOR_VERSION}.${GMP_PATCHLEVEL_VERSION})
    endif()
  endforeach()

  set(GMP_VERSION
    ${GMP_MAJOR_VERSION}.${GMP_MINOR_VERSION}.${GMP_PATCHLEVEL_VERSION})

  # Check whether found version exists and exceeds the minimum requirement
  if(NOT GMP_VERSION)
    set(GMP_VERSION_OK FALSE)
    message(STATUS "GMP version was not detected")
  elseif(${GMP_VERSION} VERSION_LESS ${GMP_FIND_VERSION})
    set(GMP_VERSION_OK FALSE)
    message(STATUS "GMP version ${GMP_VERSION} found in ${GMP_INCLUDE_DIRS}, "
                   "but at least version ${GMP_FIND_VERSION} is required")
  else()
    set(GMP_VERSION_OK TRUE)
  endif()
endmacro(_gmp_check_version)

if(NOT GMP_VERSION_OK)
  set(GMP_INCLUDE_DIRS NOTFOUND)
  set(GMP_LIBRARIES NOTFOUND)

  # search first if an GMPConfig.cmake is available in the system,
  # if successful this would set GMP_INCLUDE_DIRS and the rest of
  # the script will work as usual
  find_package(GMP ${GMP_FIND_VERSION} NO_MODULE QUIET)

  if(NOT GMP_INCLUDE_DIRS)
    find_path(GMP_INCLUDE_DIRS NAMES gmp.h
      HINTS ENV GMPDIR ENV GMPDIR
      PATHS ${INCLUDE_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/include
      )
  endif()

  if(GMP_INCLUDE_DIRS)
    _GMP_check_version()
  endif()

  if(NOT GMP_LIBRARIES)
    find_library(GMP_LIBRARIES NAMES gmp
      HINTS ENV GMPDIR ENV GMPDIR
      PATHS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
      )
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIRS GMP_LIBRARIES GMP_VERSION_OK)

  mark_as_advanced(GMP_INCLUDE_DIRS GMP_LIBRARIES)

endif()
