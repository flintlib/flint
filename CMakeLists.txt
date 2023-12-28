#
# Copyright (C) 2023 Isuru Fernando
# Copyright (C) 2023 Albin Ahlbäck
#
# This file is part of FLINT.
#
# FLINT is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License (LGPL) as published
# by the Free Software Foundation; either version 2.1 of the License, or
# (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 3.22)

include(CheckCCompilerFlag)
include(CheckCSourceRuns)
include(CheckIPOSupported)
include(FindPkgConfig)

# Source of truth for project version
file(STRINGS VERSION FLINT_VERSION_FULL)
string(REGEX MATCH "([0-9]+\.[0-9]+\.[0-9]+)" _ ${FLINT_VERSION_FULL})
set(FLINT_VERSION ${CMAKE_MATCH_1})

project(flint
  VERSION ${FLINT_VERSION}
  DESCRIPTION "Fast Library for Number Theory"
  HOMEPAGE_URL https://flintlib.org/
  LANGUAGES C CXX)

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/configure.ac" CONFIGURE_CONTENTS)
foreach(version in MAJOR MINOR PATCH)
  # Set soname version for the library
  string(REGEX MATCH "FLINT_${version}_SO=([0-9]*)" _ ${CONFIGURE_CONTENTS})
  set(FLINT_${version}_SO ${CMAKE_MATCH_1})
  # Set flint version for the header
  set(FLINT_${version} ${PROJECT_VERSION_${version}})
endforeach()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build" FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Try to enable the fft_small module
if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    check_c_compiler_flag("-march=native" HAS_FLAG_GCC_MARCH_NATIVE)
elseif("${CMAKE_C_COMPILER_ID}" MATCHES "MSVC")
    # Check if AVX2 is available
    check_c_compiler_flag("/arch:AVX2" HAS_FLAG_MSVC_AVX2)
endif()

# Build options
option(BUILD_SHARED_LIBS "Build shared libs" on)
option(WITH_NTL "Build tests for NTL interface or not" off)

# Check if strongly ordered memory
set(STRONGLY_ORDERED_CPUS x86_64 x86 i386 i586 AMD64)
if(CMAKE_SYSTEM_PROCESSOR IN_LIST STRONGLY_ORDERED_CPUS)
    message(STATUS "Checking if system is strongly ordered - yes")
    set(FLINT_KNOW_STRONG_ORDER ON)
else()
    message(STATUS "Checking if system is strongly ordered - unsure")
    set(FLINT_KNOW_STRONG_ORDER OFF)
endif()

# Find dependencies
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp>=6.2.1)
pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr>=4.1.0)
if (WITH_NTL)
    find_package(NTL REQUIRED)
endif()
find_package(Python REQUIRED)

find_package(CBLAS)
set(FLINT_USES_BLAS ${CBLAS_FOUND})

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(FLINT_WANT_ASSERT ON)
endif()

# pthread configuration

if(MSVC)
    # Prefer vcpkg's pthreads
    find_package(PThreads4W)
    if(PThreads4W_FOUND)
        set(PThreads_LIBRARIES PThreads4W::PThreads4W)
        set(PThreads_INCLUDE_DIRS ${PThreads4W_INCLUDE_DIRS})
    else()
        find_package(PThreads REQUIRED)
    endif()
    set(FLINT_USES_PTHREAD ON CACHE BOOL "Use POSIX Threads.")
else()
    option(CMAKE_THREAD_PREFER_PTHREAD "Prefer pthreads" yes)
    option(THREADS_PREFER_PTHREAD_FLAG "Prefer -pthread flag" yes)
    find_package(Threads REQUIRED)
    set(PThreads_LIBRARIES Threads::Threads)
    set(FLINT_USES_PTHREAD ON CACHE BOOL "Use POSIX Threads.")
endif()


# Check if fft_small module is available
message(STATUS "Checking whether fft_small module is available")
set(CMAKE_REQUIRED_INCLUDES ${GMP_INCLUDE_DIRS})
if(HAS_FLAG_GCC_MARCH_NATIVE)
    set(CMAKE_REQUIRED_FLAGS "-march=native")
elseif(HAS_FLAG_MSVC_AVX2)
    set(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
endif()
check_c_source_compiles([[
    #include <gmp.h>
    #if GMP_LIMB_BITS != 64
    # error
    error
    #endif

    #include <arm_neon.h>
    #if !(defined(__GNUC__) && defined(__ARM_NEON))
    # if !(defined(_MSC_VER) && defined(_M_ARM64))
    #  error
    error
    # endif
    #endif
    void main(){};]] FLINT_FFT_SMALL_ARM)
check_c_source_compiles([[
    #include <gmp.h>
    #if GMP_LIMB_BITS != 64
    # error
    error
    #endif

    #if defined(__GNUC__)
    # include <x86intrin.h>
    #elif defined(_MSC_VER)
    # include <intrin.h>
    #else
    # error
    error
    #endif

    #if !defined(__AVX2__)
    # error
    error
    #endif
    void main(){};]] FLINT_FFT_SMALL_X86)
unset(CMAKE_REQUIRED_INCLUDES)
unset(CMAKE_REQUIRED_FLAGS)

if(FLINT_FFT_SMALL_ARM OR FLINT_FFT_SMALL_X86)
    message(STATUS "Checking whether fft_small module is available - yes")
    set(FFT_SMALL fft_small)
    set(FLINT_HAVE_FFT_SMALL ON)
else()
    message(STATUS "Checking whether fft_small module is available - no")
    set(FFT_SMALL)
endif()

# Find sources

set(_BUILD_DIRS
    generic_files
    thread_pool                     thread_support

    ulong_extras
    long_extras
    perm
    double_extras   d_vec           d_mat
    mpn_extras
    mpf_vec                         mpf_mat
    mpfr_vec                        mpfr_mat
    nmod            nmod_vec        nmod_mat        nmod_poly
    fmpz            fmpz_vec        fmpz_mat        fmpz_poly
    fmpz_mod        fmpz_mod_vec    fmpz_mod_mat    fmpz_mod_poly
    fmpq            fmpq_vec        fmpq_mat        fmpq_poly

    fq              fq_vec          fq_mat          fq_poly
    fq_nmod         fq_nmod_vec     fq_nmod_mat     fq_nmod_poly
    fq_zech         fq_zech_vec     fq_zech_mat     fq_zech_poly
    fq_default                      fq_default_mat  fq_default_poly
    fq_embed
    fq_nmod_embed
    fq_zech_embed
    padic                           padic_mat       padic_poly
    qadic

    nmod_poly_factor                fmpz_factor
    fmpz_poly_factor                fmpz_mod_poly_factor
    fq_poly_factor                  fq_nmod_poly_factor
    fq_zech_poly_factor             fq_default_poly_factor

    nmod_poly_mat                   fmpz_poly_mat

    mpoly           nmod_mpoly      fmpz_mpoly      fmpz_mod_mpoly
    fmpq_mpoly      fq_nmod_mpoly   fq_zech_mpoly

    nmod_mpoly_factor               fmpz_mpoly_factor
    fmpz_mod_mpoly_factor           fmpq_mpoly_factor
    fq_nmod_mpoly_factor            fq_zech_mpoly_factor

    fft             ${FFT_SMALL}    fmpz_poly_q    fmpz_lll
    n_poly          arith           qsieve          aprcl

    nf              nf_elem         qfb

    double_interval dlog
    fmpz_extras     fmpzi
    bool_mat        partitions
    mag
    arf             acf             arb             acb
    arb_mat         arb_poly        arb_calc        arb_hypgeom
    acb_mat         acb_poly        acb_calc        acb_hypgeom
    arb_fmpz_poly   arb_fpwrap
    acb_dft         acb_elliptic    acb_modular     acb_dirichlet
    acb_theta       dirichlet       bernoulli       hypgeom

    gr              gr_generic      gr_vec          gr_mat
    gr_poly         gr_mpoly        gr_special

    calcium
    fmpz_mpoly_q
    fexpr           fexpr_builtin
    qqbar
    ca              ca_ext          ca_field        ca_vec
    ca_poly         ca_mat

)
string(REGEX REPLACE "([A-Za-z0-9_-]+;|[A-Za-z0-9_-]+$)" "src/\\1" BUILD_DIRS "${_BUILD_DIRS}")

# NOTE: Template directories are not supposed to be used.

set(_HEADERS
    NTL-interface.h flint.h longlong.h flint-config.h gmpcompat.h fft_tuning.h
    profiler.h templates.h exception.h
)
string(REGEX REPLACE "([A-Za-z0-9_-]+\.h;|[A-Za-z0-9_-]+\.h$)" "src/\\1" HEADERS "${_HEADERS}")

# Setup for flint-config.h
check_c_compiler_flag("-mpopcnt" HAS_FLAG_MPOPCNT)
check_c_compiler_flag("-funroll-loops" HAS_FLAG_UNROLL_LOOPS)

if(HAS_FLAG_MPOPCNT)
  set(CMAKE_REQUIRED_FLAGS "-mpopcnt")
endif()
# Testing __builtin_popcountl...
check_c_source_runs([[int main(int argc, char ** argv) {
  #if defined(_WIN64)
  return __builtin_popcountll(argc) == 100;
  #else
  return __builtin_popcountl(argc) == 100;
  #endif
  }]] FLINT_USES_POPCNT)
unset(CMAKE_REQUIRED_FLAGS)

# fenv configuration
check_c_source_compiles([[#include <fenv.h>
  #ifndef FE_DOWNWARD
  # error FE_DOWNWARD not available
  #endif
  void main(){};]] FLINT_USES_FENV)

# cpu_set_t configuration
set(CMAKE_REQUIRED_FLAGS "${PThreads_LIBRARIES}")
check_c_source_compiles([[#define _GNU_SOURCE
  #include <sched.h>
  #include <pthread.h>
  int main() { cpu_set_t s; CPU_ZERO(&s);
  pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), 0);
  return 0; }]] FLINT_USES_CPUSET)
unset(CMAKE_REQUIRED_FLAGS)

# Thread-local storage configuration
set(FLINT_USES_TLS ON CACHE BOOL "Use thread local storage.")

# Memory manager configuration
set(MEMORY_MANAGER "single" CACHE STRING "The FLINT memory manager.")
set_property(CACHE MEMORY_MANAGER PROPERTY STRINGS single reentrant gc)
message(STATUS "Using FLINT memory manager: ${MEMORY_MANAGER}")

if(MEMORY_MANAGER STREQUAL "reentrant")
	set(FLINT_REENTRANT ON)
else()
	set(FLINT_REENTRANT OFF)
endif()

# gmpcompat.h configuration
set(CMAKE_REQUIRED_INCLUDES ${GMP_INCLUDE_DIRS})
check_c_source_compiles([[#include <gmp.h>
  #ifndef _LONG_LONG_LIMB
  # error mp_limb_t != unsigned long long limb
  #endif
  void main(){};]] LONG_LONG_LIMB)

# Populate headers
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/CMake/cmake_config.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/src/flint-config.h
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/flint.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/src/flint.h
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fmpz/link/fmpz_${MEMORY_MANAGER}.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fmpz/fmpz.c
    COPYONLY
)
if(LONG_LONG_LIMB)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gmpcompat-longlong.h.in
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gmpcompat.h
        COPYONLY
    )
else()
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gmpcompat.h.in
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gmpcompat.h
        COPYONLY
    )
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/src/fft_tuning64.in
        ${CMAKE_CURRENT_SOURCE_DIR}/src/fft_tuning.h
        COPYONLY
    )   
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/src/fft_tuning32.in
        ${CMAKE_CURRENT_SOURCE_DIR}/src/fft_tuning.h
        COPYONLY
    )
endif()


foreach (build_dir IN LISTS BUILD_DIRS)
    file(GLOB TEMP RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${build_dir}/*.c")
    list(APPEND SOURCES ${TEMP})
    file(GLOB TEMP RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${build_dir}/*.h")
    list(APPEND HEADERS ${TEMP})
endforeach ()

set(TEMP ${HEADERS})
set(HEADERS )
foreach(header IN LISTS TEMP)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${header})
        list(APPEND HEADERS ${header})
    else()
        list(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${header})  
    endif()
endforeach()

file(GLOB TEMP "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
list(APPEND HEADERS ${TEMP})

add_library(flint ${SOURCES})
target_link_libraries(flint PUBLIC
    PkgConfig::GMP PkgConfig::MPFR ${PThreads_LIBRARIES})

if(FLINT_USES_BLAS)
    target_link_libraries(flint PUBLIC ${CBLAS_LIBRARIES})
endif()

# Include directories

target_include_directories(flint PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src ${GMP_INCLUDE_DIRS} ${MPFR_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR} ${PThreads_INCLUDE_DIRS})

if(FLINT_USES_BLAS)
    target_include_directories(flint PUBLIC ${CBLAS_INCLUDE_DIRS})
endif()

if(BUILD_SHARED_LIBS AND MSVC)
    # Export all functions automatically (except global data)
    set_target_properties(flint PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    # Export flint's global data that are marked manually
    target_compile_definitions(flint PRIVATE "FLINT_BUILD_DLL")
endif()

if (HAS_FLAG_MPOPCNT)
    target_compile_options(flint PUBLIC "-mpopcnt")
endif()
if (HAS_FLAG_UNROLL_LOOPS)
    target_compile_options(flint PUBLIC "-funroll-loops")
endif()

if(HAS_FLAG_GCC_MARCH_NATIVE)
    target_compile_options(flint PUBLIC "-march=native")
elseif(HAS_FLAG_MSVC_AVX2)
    target_compile_options(flint PUBLIC "/arch:AVX2")
endif()

# Versioning

set_target_properties(flint PROPERTIES
    VERSION ${FLINT_MAJOR_SO}.${FLINT_MINOR_SO}.${FLINT_PATCH_SO}
    SOVERSION ${FLINT_MAJOR_SO}
)

# Following versioning parts are optional
# Match versioning scheme in configure based build system.
if (APPLE)
    if(${CMAKE_VERSION} VERSION_LESS "3.17.0")
        message(WARNING "To match the versioning scheme of configure based build system, switch to cmake 3.17.0")
    else ()
        set_target_properties(flint PROPERTIES
            MACHO_COMPATIBILITY_VERSION ${FLINT_MAJOR_SO}.${FLINT_MINOR_SO}
            MACHO_CURRENT_VERSION ${FLINT_MAJOR_SO}.${FLINT_MINOR_SO}.${FLINT_PATCH_SO}
        )
    endif()
elseif (WIN32)
    set_target_properties(flint PROPERTIES RUNTIME_OUTPUT_NAME "flint-${FLINT_MAJOR_SO}")
endif()

set_property(TARGET flint PROPERTY C_STANDARD 11)

if(NOT DEFINED IPO_SUPPORTED)
    message(STATUS "Checking for IPO")
    check_ipo_supported(RESULT ipo_supported LANGUAGES C)
    if(ipo_supported)
	if (MSVC)
            message(STATUS "Checking for IPO - found, but disabled for MSVC")
            set(ipo_supported FALSE)
        else()
            message(STATUS "Checking for IPO - found")
        endif()
    else()
        message(STATUS "Checking for IPO - not found")
    endif()
    set(IPO_SUPPORTED ${ipo_supported} CACHE INTERNAL "Introprocedural Optimization" FORCE)
endif()

# allow overriding IPO by setting -DCMAKE_INTERPROCEDURAL_OPTIMIZATION
if (IPO_SUPPORTED AND "${CMAKE_INTERPROCEDURAL_OPTIMIZATION}" STREQUAL "")
    set_target_properties(flint PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

if(NOT MSVC)
	target_link_libraries(flint PUBLIC m)
endif()

include(GNUInstallDirs)

install(TARGETS flint
            RUNTIME DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}"
            ARCHIVE DESTINATION "${CMAKE_INSTALL_FULL_LIBDIR}"
            LIBRARY DESTINATION "${CMAKE_INSTALL_FULL_LIBDIR}"
        )

install(FILES ${HEADERS} DESTINATION include/flint)

set_target_properties(flint
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)

if(BUILD_TESTING)
    set(FLINT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src) # To get src/test/main
    enable_testing()
    foreach(dir IN LISTS BUILD_DIRS FLINT_SRC)
        file(GLOB test_file "${dir}/test/main.c")
        if(test_file)
            file(RELATIVE_PATH test_name ${CMAKE_CURRENT_SOURCE_DIR} ${test_file})
            string(REPLACE "/" "-" test_name ${test_name})
            get_filename_component(test_name ${test_name} NAME_WE)
            add_executable(${test_name} ${test_file})
            target_link_libraries(${test_name} flint)

            add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)

            set_target_properties(${test_name}
                PROPERTIES
                ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
        endif()
    endforeach ()

    if(WITH_NTL)
        file(GLOB test_file "src/interfaces/test/t-NTL-interface.cpp")
        file(RELATIVE_PATH test_name ${CMAKE_CURRENT_SOURCE_DIR} ${test_file})
        string(REPLACE "/" "-" test_name ${test_name})
        get_filename_component(test_name ${test_name} NAME_WE)
        add_executable(${test_name} ${test_file})
        target_link_libraries(${test_name} flint ${NTL_LIBRARY})
        target_include_directories(${test_name} PUBLIC ${NTL_INCLUDE_DIR})

        add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)

        set_target_properties(${test_name}
            PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
    endif()
endif()


if(BUILD_DOCS)
    find_package(Sphinx REQUIRED)
    file(GLOB DOC_SOURCES doc/source/*.rst)
    add_custom_target(html
        COMMAND ${SPHINX_EXECUTABLE} -b html "${CMAKE_CURRENT_SOURCE_DIR}/doc/source" "${CMAKE_CURRENT_BINARY_DIR}/html"
        SOURCES ${DOC_SOURCES})  
    add_custom_target(latex
        COMMAND ${SPHINX_EXECUTABLE} -b latex "${CMAKE_CURRENT_SOURCE_DIR}/doc/source" "${CMAKE_CURRENT_BINARY_DIR}/latex"
        SOURCES ${DOC_SOURCES})  
    add_custom_target(pdf DEPENDS latex COMMAND make WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/latex") 
endif()
