
set(GMP_ROOT_DIR "${GMP_ROOT_DIR}"  CACHE PATH "Directory to search for gmp" )

find_package(PkgConfig QUIET)
if( PkgConfig_FOUND )
  pkg_search_module(PC_GMP QUIET gmp)
  if( PC_GMP_FOUND )
    set( GMP_VERSION ${PC_GMP_VERSION} )
  endif()
endif()

find_path( GMP_INCLUDE_DIR
  NAMES gmp.h
  PATHS "${GMP_ROOT_DIR}"
  HINTS ${PC_GMP_INCLUDEDIR} ${PC_GMP_INCLUDE_DIRS}
  )
find_library( GMP_LIBRARY
  NAMES gmp
  PATHS "${GMP_ROOT_DIR}"
  HINTS ${PC_GMP_LIBDIR} ${PC_GMP_LIBRARY_DIRS}
  )

if(NOT PC_GMP_FOUND)
  set( _VERSION_FILE ${GMP_INCLUDE_DIR}/gmp.h )
  if( EXISTS ${_VERSION_FILE} )
    file( STRINGS ${_VERSION_FILE} _VERSION_LINE REGEX "define[ ]+__GNU_MP_VERSION[ ]+" )
    if( _VERSION_LINE )
      string( REGEX REPLACE ".*define[ ]+__GNU_MP_VERSION[ ]+(.*)[ ]*" "\\1" GMP_VERSION_MAJOR "${_VERSION_LINE}" )
    endif()
    file( STRINGS ${_VERSION_FILE} _VERSION_LINE REGEX "define[ ]+__GNU_MP_VERSION_MINOR[ ]+" )
    if( _VERSION_LINE )
      string( REGEX REPLACE ".*define[ ]+__GNU_MP_VERSION_MINOR[ ]+(.*)[ ]*" "\\1" GMP_VERSION_MINOR "${_VERSION_LINE}" )
    endif()
    file( STRINGS ${_VERSION_FILE} _VERSION_LINE REGEX "define[ ]+__GNU_MP_VERSION_PATCHLEVEL[ ]+" )
    if( _VERSION_LINE )
      string( REGEX REPLACE ".*define[ ]+__GNU_MP_VERSION_PATCHLEVEL[ ]+(.*)[ ]*" "\\1" GMP_VERSION_PATCHLEVEL "${_VERSION_LINE}" )
    endif()
    set( GMP_VERSION "${GMP_VERSION_MAJOR}.${GMP_VERSION_MINOR}.${GMP_VERSION_PATCHLEVEL}")
  endif()
  unset( _VERSION_FILE )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( gmp
  FOUND_VAR GMP_FOUND
  REQUIRED_VARS
    GMP_LIBRARY
    GMP_INCLUDE_DIR
  VERSION_VAR GMP_VERSION
  )

if(GMP_FOUND)
  set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
  set(GMP_LIBRARIES ${GMP_LIBRARY})
  if(NOT TARGET gmp::gmp)
    add_library(gmp::gmp UNKNOWN IMPORTED)
    set_target_properties( gmp::gmp
      PROPERTIES
        IMPORTED_LOCATION "${GMP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}" 
      )
  endif()
  mark_as_advanced(GMP_ROOT_DIR)
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
