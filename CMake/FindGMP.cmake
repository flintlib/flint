# Try to find the GMP library
# https://gmplib.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(GMP 6.2.1)
# to require version 6.2.1 to newer of GMP.
#
# Once done this will define
#
#  GMP_FOUND - system has GMP lib with correct version
#  GMP_INCLUDE_DIRS - the GMP include directory
#  GMP_LIBRARIES - the GMP library
#

# Set GMP_FIND_VERSION to 6.0.0 if no minimum version is specified
if(NOT GMP_FIND_VERSION)
  if(NOT GMP_FIND_VERSION_MAJOR)
    set(GMP_FIND_VERSION_MAJOR 6)
  endif()
  if(NOT GMP_FIND_VERSION_MINOR)
    set(GMP_FIND_VERSION_MINOR 0)
  endif()
  if(NOT GMP_FIND_VERSION_PATCH)
    set(GMP_FIND_VERSION_PATCH 0)
  endif()
  set(GMP_FIND_VERSION
    "${GMP_FIND_VERSION_MAJOR}.${GMP_FIND_VERSION_MINOR}.${GMP_FIND_VERSION_PATCH}")
endif()

find_path(GMP_INCLUDE_DIRS
          NAMES gmp.h
          PATHS $ENV{GMPDIR} ${INCLUDE_INSTALL_DIR})

find_library(GMP_LIBRARIES
             gmp
             PATHS $ENV{GMPDIR} ${LIB_INSTALL_DIR})

if(GMP_INCLUDE_DIRS AND GMP_LIBRARIES)

  # This program will fail to compile if GMP is too old.
  # We prefer to perform this "test" at compile-time to
  # avoid problems with e.g. try_run() during cross-compilation.
  file(WRITE ${PROJECT_BINARY_DIR}/gmp-version-check.c ""
  "#include <gmp.h>\n"
  "\n"
  "#define GMP_FIND_VERSION_MAJOR ${GMP_FIND_VERSION_MAJOR}\n"
  "#define GMP_FIND_VERSION_MINOR ${GMP_FIND_VERSION_MINOR}\n"
  "#define GMP_FIND_VERSION_PATCH ${GMP_FIND_VERSION_PATCH}\n"
  "\n"
  "#if __GNU_MP_VERSION < GMP_FIND_VERSION_MAJOR\n"
  "#error insufficient GMP major version\n"
  "#elif __GNU_MP_VERSION == GMP_FIND_VERSION_MAJOR\n"
  "#if __GNU_MP_VERSION_MINOR < GMP_FIND_VERSION_MINOR\n"
  "#error insufficient GMP minor version\n"
  "#elif __GNU_MP_VERSION_MINOR == GMP_FIND_VERSION_MINOR\n"
  "#if __GNU_MP_VERSION_PATCH < GMP_FIND_VERSION_PATCH\n"
  "#error insufficient GMP patch version\n"
  "#endif\n"
  "#endif\n"
  "#endif\n"
  "\n"
  "int main(int argc, char** argv) { return 0; }\n")

  # Try to compile the test program above with the appropriate version
  # strings substituted in.
  try_compile(GMP_VERSION_OK
          "${PROJECT_BINARY_DIR}"
          "${PROJECT_BINARY_DIR}/gmp-version-check.c"
          CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${GMP_INCLUDE_DIRS}")
endif()

if(NOT GMP_VERSION_OK)
  message(STATUS "No sufficient GMP version detected")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG
                                  GMP_INCLUDE_DIRS GMP_LIBRARIES GMP_VERSION_OK)
mark_as_advanced(GMP_INCLUDE_DIRS GMP_LIBRARIES)
