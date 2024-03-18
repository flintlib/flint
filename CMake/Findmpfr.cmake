
set(MPFR_ROOT_DIR "${MPFR_ROOT_DIR}"  CACHE PATH "Directory to search for mpfr" )

find_package(PkgConfig QUIET)
if( PkgConfig_FOUND )
  pkg_search_module(PC_MPFR QUIET mpfr)
  if( PC_MPFR_FOUND )
    set( MPFR_VERSION ${PC_MPFR_VERSION} )
  endif()
endif()

find_path( MPFR_INCLUDE_DIR
  NAMES mpfr.h
  PATHS "${MPFR_ROOT_DIR}"
  HINTS ${PC_MPFR_INCLUDEDIR} ${PC_MPFR_INCLUDE_DIRS}
  )
find_library( MPFR_LIBRARY
  NAMES mpfr
  PATHS "${MPFR_ROOT_DIR}"
  HINTS ${PC_MPFR_LIBDIR} ${PC_MPFR_LIBRARY_DIRS}
  )

if(NOT PC_MPFR_FOUND)
  set( _VERSION_FILE ${MPFR_INCLUDE_DIR}/mpfr.h )
  if( EXISTS ${_VERSION_FILE} )
    file( STRINGS ${_VERSION_FILE} _VERSION_LINE REGEX "define[ ]+MPFR_VERSION_STRING" )
    if( _VERSION_LINE )
      string( REGEX REPLACE ".*define[ ]+MPFR_VERSION_STRING[ ]+\"([^\"]*)\".*" "\\1" MPFR_VERSION "${_VERSION_LINE}" )
    endif()
  endif()
  unset( _VERSION_FILE )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( mpfr
  FOUND_VAR MPFR_FOUND
  REQUIRED_VARS
    MPFR_LIBRARY
    MPFR_INCLUDE_DIR
  VERSION_VAR MPFR_VERSION
  )

if(MPFR_FOUND)
  set(MPFR_INCLUDE_DIRS ${MPFR_INCLUDE_DIR})
  set(MPFR_LIBRARIES ${MPFR_LIBRARY})
  if(NOT TARGET mpfr::mpfr)
    add_library(mpfr::mpfr UNKNOWN IMPORTED)
    set_target_properties( mpfr::mpfr
      PROPERTIES
        IMPORTED_LOCATION "${MPFR_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}" 
      )
  endif()
  mark_as_advanced(MPFR_ROOT_DIR)
endif()

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)
