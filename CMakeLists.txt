cmake_minimum_required(VERSION 3.9)

include(CheckIPOSupported)
include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckPrototypeDefinition)
include(CheckCCompilerFlag)
include(CheckCSourceRuns)
include(CheckLibraryExists)
include(TestBigEndian)

project(flint LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build" FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(BUILD_SHARED_LIBS "Build shared libs" on)
option(WITH_NTL "Build with NTL or not" off)

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/configure" CONFIGURE_CONTENTS)
string(REGEX MATCH "FLINT_MAJOR=([0-9]*)" _ ${CONFIGURE_CONTENTS})
set(FLINT_MAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "FLINT_MINOR=([0-9]*)" _ ${CONFIGURE_CONTENTS})
set(FLINT_MINOR ${CMAKE_MATCH_1})
string(REGEX MATCH "FLINT_PATCH=([0-9]*)" _ ${CONFIGURE_CONTENTS})
set(FLINT_PATCH ${CMAKE_MATCH_1})

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
if (WITH_NTL)
    find_package(NTL REQUIRED)
endif()
find_package(PythonInterp REQUIRED)

find_package(CBLAS)
set(FLINT_USES_BLAS ${CBLAS_FOUND})

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(FLINT_WANT_ASSERT ON)
endif()

# pthread configuration

if(MSVC)
    find_package(PThreads REQUIRED)
    set(FLINT_USES_PTHREAD ON CACHE BOOL "Use POSIX Threads.")
else()
    option(CMAKE_THREAD_PREFER_PTHREAD "Prefer pthreads" yes)
    option(THREADS_PREFER_PTHREAD_FLAG "Prefer -pthread flag" yes)
    find_package(Threads REQUIRED)
    set(PThreads_LIBRARIES Threads::Threads)
    set(FLINT_USES_PTHREAD ON CACHE BOOL "Use POSIX Threads.")
endif()

# Find sources

set(BUILD_DIRS
    aprcl ulong_extras long_extras perm fmpz fmpz_vec fmpz_poly 
    fmpq_poly fmpz_mat fmpz_lll mpfr_vec mpfr_mat mpf_vec mpf_mat nmod_vec nmod_poly 
    nmod_poly_factor arith mpn_extras nmod_mat fmpq fmpq_vec fmpq_mat padic 
    fmpz_poly_q fmpz_poly_mat nmod_poly_mat fmpz_mod_poly fmpz_mod_mat 
    fmpz_mod_poly_factor fmpz_factor fmpz_poly_factor fft qsieve 
    double_extras d_vec d_mat padic_poly padic_mat qadic  
    fq fq_vec fq_mat fq_poly fq_poly_factor
    fq_nmod fq_nmod_vec fq_nmod_mat fq_nmod_poly fq_nmod_poly_factor 
    fq_zech fq_zech_vec fq_zech_mat fq_zech_poly fq_zech_poly_factor 
    fq_default fq_default_mat fq_default_poly fq_default_poly_factor
    thread_pool thread_support fmpz_mod fmpz_mod_vec n_poly
    mpoly fmpz_mpoly fmpq_mpoly nmod_mpoly fq_nmod_mpoly fmpz_mod_mpoly 
    fmpz_mpoly_factor fmpq_mpoly_factor nmod_mpoly_factor fmpz_mod_mpoly_factor 
    fq_nmod_mpoly_factor fq_zech_mpoly fq_zech_mpoly_factor 
    fq_embed fq_nmod_embed fq_zech_embed nmod
)

set(TEMPLATE_DIRS
    fq_vec_templates fq_mat_templates fq_poly_templates
    fq_poly_factor_templates fq_embed_templates fq_templates
)

set(SOURCES
    printf.c fprintf.c sprintf.c scanf.c fscanf.c sscanf.c clz_tab.c
    memory_manager.c version.c profiler.c exception.c
    hashmap.c inlines.c fmpz/fmpz.c
)

if (MSVC)
    list(APPEND SOURCES gettimeofday.c)
endif()

if (WITH_NTL)
    list(APPEND SOURCES interfaces/NTL-interface.cpp)
endif()

set(HEADERS
    NTL-interface.h flint.h longlong.h flint-config.h gmpcompat.h fft_tuning.h
    fmpz-conversions.h profiler.h templates.h exception.h hashmap.h
)

foreach (build_dir IN LISTS BUILD_DIRS TEMPLATE_DIRS)
    file(GLOB TEMP RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${build_dir}/*.c")
    list(APPEND SOURCES ${TEMP})
    file(GLOB TEMP RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${build_dir}.h")
    list(APPEND HEADERS ${TEMP})
endforeach ()

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c 
"
from os.path import join

with open(join('${CMAKE_CURRENT_SOURCE_DIR}','qadic', 'CPimport.txt')) as fin:
    with open(join('${CMAKE_CURRENT_BINARY_DIR}','CPimport.h.in'), 'w+') as fout:
        while True:
            l = fin.readline()
            if not l:
                break
            l = l.replace(' ', ',')
            l = l.replace('\\n', ',\\n')
            fout.writelines([l])
"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
configure_file(${CMAKE_CURRENT_BINARY_DIR}/CPimport.h.in ${CMAKE_CURRENT_BINARY_DIR}/CPimport.h COPYONLY)

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
if(cxx_thread_local IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  set(FLINT_USES_TLS ON CACHE BOOL "Use thread local storage.")
endif()

# Memory manager configuration
set(MEMORY_MANAGER "reentrant" CACHE STRING "The FLINT memory manager.")
set_property(CACHE MEMORY_MANAGER PROPERTY STRINGS single reentrant gc)
message(STATUS "Using FLINT memory manager: ${MEMORY_MANAGER}")

if(MEMORY_MANAGER STREQUAL "reentrant")
	set(FLINT_REENTRANT ON)
else()
	set(FLINT_REENTRANT OFF)
endif()

# Populate headers
configure_file(
    config.h.in
    flint-config.h
)

configure_file(
    fmpz-conversions-${MEMORY_MANAGER}.in
    fmpz-conversions.h
    COPYONLY
)

configure_file(
    fmpz/link/fmpz_${MEMORY_MANAGER}.c
    fmpz/fmpz.c
    COPYONLY
)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    configure_file(
        fft_tuning64.in
        fft_tuning.h
        COPYONLY
    )   
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
    configure_file(
        fft_tuning32.in
        fft_tuning.h
        COPYONLY
    )
endif()


set(TEMP ${HEADERS})
set(HEADERS )
foreach(header IN LISTS TEMP)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${header})
        list(APPEND HEADERS ${header})
    else()
        list(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${header})  
    endif()
endforeach()

add_library(flint ${SOURCES})
target_link_libraries(flint PUBLIC
    ${NTL_LIBRARY} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} ${PThreads_LIBRARIES}
)

if(FLINT_USES_BLAS)
    target_link_libraries(flint PUBLIC ${CBLAS_LIBRARIES})
endif()

# Include directories

target_include_directories(flint PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR} ${GMP_INCLUDE_DIRS} ${MPFR_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR} ${PThreads_INCLUDE_DIRS}
    ${NTL_INCLUDE_DIR}
)

if(FLINT_USES_BLAS)
    target_include_directories(flint PUBLIC ${CBLAS_INCLUDE_DIRS})
endif()

if(BUILD_SHARED_LIBS AND WIN32)
    # Export flint's functions
    target_compile_definitions(flint PRIVATE "FLINT_BUILD_DLL")
    # Use MPIR's dll import functions
    target_compile_definitions(flint PUBLIC "MSC_USE_DLL")
endif()

if (HAS_FLAG_MPOPCNT)
    target_compile_options(flint PUBLIC "-mpopcnt")
endif()
if (HAS_FLAG_UNROLL_LOOPS)
    target_compile_options(flint PUBLIC "-funroll-loops")
endif()

# Versioning

set_target_properties(flint PROPERTIES
    VERSION ${FLINT_MAJOR}.${FLINT_MINOR}.${FLINT_PATCH}
    SOVERSION ${FLINT_MAJOR}
)

# Following versioning parts are optional
# Match versioning scheme in configure based build system.
if (APPLE)
    if(${CMAKE_VERSION} VERSION_LESS "3.17.0")
        message(WARNING "To match the versioning scheme of configure based build system, switch to cmake 3.17.0")
    else ()
        set_target_properties(flint PROPERTIES
            MACHO_COMPATIBILITY_VERSION ${FLINT_MAJOR}.${FLINT_MINOR}
            MACHO_CURRENT_VERSION ${FLINT_MAJOR}.${FLINT_MINOR}.${FLINT_PATCH}
        )
    endif()
elseif (WIN32)
    set_target_properties(flint PROPERTIES RUNTIME_OUTPUT_NAME "flint-${FLINT_MAJOR}")
endif()

if(NOT DEFINED IPO_SUPPORTED)
    message(STATUS "Checking for IPO")
    check_ipo_supported(RESULT ipo_supported LANGUAGES C)
    if(ipo_supported)
        message(STATUS "Checking for IPO - found")
    else()
        message(STATUS "Checking for IPO - not found")
    endif()
    set(IPO_SUPPORTED ${ipo_supported} CACHE INTERNAL "Introprocedural Optimization" FORCE)
endif()

if(IPO_SUPPORTED)
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

file(GLOB TEMP "flintxx/*.h")
install(FILES ${TEMP} DESTINATION include/flint/flintxx)

file(GLOB TEMP "*xx.h")
install(FILES ${TEMP} DESTINATION include/flint)

set_target_properties(flint
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)

if(BUILD_TESTING)
    enable_testing()
    add_library(test_helpers STATIC test_helpers.c)
    target_link_libraries(test_helpers flint)

    foreach (build_dir IN LISTS BUILD_DIRS CMAKE_CURRENT_SOURCE_DIR)
        file(GLOB TEST_FILES "${build_dir}/test/*.c")
        foreach(test_file IN LISTS TEST_FILES)
            file(RELATIVE_PATH test_name ${CMAKE_CURRENT_SOURCE_DIR} ${test_file})
            string(REPLACE "/" "-" test_name ${test_name})
            get_filename_component(test_name ${test_name} NAME_WE)
            add_executable(${test_name} ${test_file})
            target_link_libraries(${test_name}
                flint test_helpers
            )

            add_test(
                NAME ${test_name}
                COMMAND $<TARGET_FILE:${test_name}>
            )

            set_target_properties(${test_name}
                PROPERTIES
                ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
            )
        endforeach()
    endforeach ()
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
