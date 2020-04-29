cmake_minimum_required(VERSION 3.9)

include(CheckIPOSupported)
include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckPrototypeDefinition)
include(CheckCCompilerFlag)
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

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
if (WITH_NTL)
    find_package(NTL REQUIRED)
endif()
find_package(PythonInterp REQUIRED)


if(MSVC)
    find_package(PThreads REQUIRED)
else()
    option(CMAKE_THREAD_PREFER_PTHREAD "Prefer pthreads" yes)
    option(THREADS_PREFER_PTHREAD_FLAG "Prefer -pthread flag" yes)
    find_package(Threads REQUIRED)
    set(PThreads_LIBRARIES Threads::Threads)
endif()

# Find sources

set(BUILD_DIRS
    aprcl ulong_extras long_extras perm fmpz fmpz_vec fmpz_sparse_vec fmpz_poly 
    fmpq_poly fmpz_mat fmpz_sparse_mat fmpz_lll mpfr_vec mpfr_mat mpf_vec mpf_mat nmod_vec nmod_sparse_vec nmod_poly 
    nmod_poly_factor arith mpn_extras nmod_mat nmod_sparse_mat fmpq fmpq_vec fmpq_mat padic 
    fmpz_poly_q fmpz_poly_mat nmod_poly_mat fmpz_mod_poly fmpz_mod_mat
    fmpz_mod_poly_factor fmpz_factor fmpz_poly_factor fft qsieve 
    double_extras d_vec d_mat padic_poly padic_mat qadic  
    fq fq_vec fq_sparse_vec fq_mat fq_sparse_mat fq_poly fq_poly_factor
    fq_nmod fq_nmod_vec fq_nmod_sparse_vec fq_nmod_mat fq_nmod_sparse_mat fq_nmod_poly fq_nmod_poly_factor 
    fq_zech fq_zech_vec fq_zech_sparse_vec fq_zech_mat fq_zech_sparse_mat fq_zech_poly fq_zech_poly_factor 
    thread_pool fmpz_mod 
    mpoly fmpz_mpoly fmpq_mpoly nmod_mpoly fq_nmod_mpoly fmpz_mod_mpoly 
    flintxx
)

set(TEMPLATE_DIRS
    fq_vec_templates fq_mat_templates fq_poly_templates
    fq_poly_factor_templates fq_templates
)

set(SOURCES
    printf.c fprintf.c sprintf.c scanf.c fscanf.c sscanf.c clz_tab.c
    memory_manager.c version.c profiler.c thread_support.c exception.c
    heap.c hashmap.c inlines.c fmpz/fmpz.c
)

if (WITH_NTL)
    list(APPEND SOURCES interfaces/NTL-interface.cpp)
endif()

set(HEADERS
    NTL-interface.h flint.h longlong.h config.h gmpcompat.h fft_tuning.h
    fmpz-conversions.h profiler.h templates.h exception.h heap.h hashmap.h
)

foreach (build_dir IN LISTS BUILD_DIRS TEMPLATE_DIRS)
    file(GLOB TEMP RELATIVE "${CMAKE_SOURCE_DIR}" "${build_dir}/*.c")
    list(APPEND SOURCES ${TEMP})
    file(GLOB TEMP RELATIVE "${CMAKE_SOURCE_DIR}" "${build_dir}/*.h")
    list(APPEND HEADERS ${TEMP})
endforeach ()

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c 
"
from os.path import join

with open(join('${CMAKE_SOURCE_DIR}','qadic', 'CPimport.txt')) as fin:
    with open('CPimport.h.in', 'w+') as fout:
        while True:
            l = fin.readline()
            if not l:
                break
            l = l.replace(' ', ',')
            l = l.replace('\\n', ',\\n')
            fout.writelines([l])
"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
configure_file(${CMAKE_BINARY_DIR}/CPimport.h.in ${CMAKE_BINARY_DIR}/CPimport.h COPYONLY)

# Configuration checks
test_big_endian(HAVE_BIG_ENDIAN)

if(NOT HAVE_BIG_ENDIAN)
    set(HAVE_LITTLE_ENDIAN ON)
    set(HAVE_DOUBLE_IEEE_LITTLE_ENDIAN ON)
else()
    set(HAVE_DOUBLE_IEEE_BIG_ENDIAN ON)
endif()

# setup for config.h
check_include_files(alloca.h HAVE_ALLOCA_H)
check_include_files(dlfcn.h HAVE_DLFCN_H)
check_include_files(fcntl.h HAVE_FCNTL_H)
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(locale.h HAVE_LOCALE_H)
check_include_files(memory.h HAVE_MEMORY_H)
check_include_files(sstream HAVE_SSTREAM)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(string.h HAVE_STRING_H)
check_include_files(sys/mman.h HAVE_SYS_MMAN_H)
check_include_files(sys/param.h HAVE_SYS_PARAM_H)
check_include_files(sys/processor.h HAVE_SYS_PROCESSOR_H)
check_include_files(sys/resource.h HAVE_SYS_RESOURCE_H)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/sysctl.h HAVE_SYS_SYSCTL_H)
check_include_files(sys/syssgi.h HAVE_SYS_SYSSGI_H)
check_include_files(sys/systemcfg.h HAVE_SYS_SYSTEMCFG_H)
check_include_files(sys/times.h HAVE_SYS_TIMES_H)
check_include_files(sys/time.h HAVE_SYS_TIME_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(dlfcn.h HAVE_DLFCN_H)
check_include_files(fpu_control.h HAVE_FPU_CONTROL_H)


check_function_exists(alloca HAVE_ALLOCA)
check_function_exists(alarm HAVE_ALARM)
check_function_exists(clock HAVE_CLOCK)
check_function_exists(clock_gettime HAVE_CLOCK_GETTIME)
check_function_exists(cputime HAVE_CPUTIME)
check_function_exists(fgetc HAVE_DECL_FGETC)
check_function_exists(fscanf HAVE_DECL_FSCANF)
check_function_exists(optarg HAVE_DECL_OPTARG)
check_function_exists(ungetc HAVE_DECL_UNGETC)
check_function_exists(vfprintf HAVE_DECL_VFPRINTF)
check_function_exists(getpagesize HAVE_GETPAGESIZE)
check_function_exists(getrusage HAVE_GETRUSAGE)
check_function_exists(gettimeofday HAVE_GETTIMEOFDAY)

check_c_compiler_flag("-mpopcnt" HAS_FLAG_MPOPCNT)
check_c_compiler_flag("-funroll-loops" HAS_FLAG_UNROLL_LOOPS)

check_library_exists(${PThreads_LIBRARIES} pthread_getaffinity_np "" HAVE_CPU_SET_T)

set(HAVE_TLS ON CACHE BOOL "Use thread local storage.")

set(MEMORY_MANAGER "reentrant" CACHE STRING "The FLINT memory manager.")
set_property(CACHE MEMORY_MANAGER PROPERTY STRINGS single reentrant gc)

if(MEMORY_MANAGER STREQUAL "reentrant")
	set(FLINT_REENTRANT ON)
else()
	set(FLINT_REENTRANT OFF)
endif()

configure_file(
    config.h.in
    config.h
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
    if(EXISTS ${CMAKE_SOURCE_DIR}/${header})
        list(APPEND HEADERS ${header})
    else()
        list(APPEND HEADERS ${CMAKE_BINARY_DIR}/${header})  
    endif()
endforeach()

file(GLOB TEMP "${CMAKE_SOURCE_DIR}/*.h")
list(APPEND HEADERS ${TEMP})

add_library(flint ${SOURCES})
target_link_libraries(flint PUBLIC
    ${NTL_LIBRARY} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} ${PThreads_LIBRARIES}
)

# Include directories

target_include_directories(flint PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR} ${GMP_INCLUDE_DIRS} ${MPFR_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR} ${PThreads_INCLUDE_DIRS}
    ${NTL_INCLUDE_DIR}
)


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

install(TARGETS flint
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
        )

install(FILES ${HEADERS} DESTINATION include/flint)

set_target_properties(flint
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

if(BUILD_TESTING)
    enable_testing()
    add_library(test_helpers STATIC test_helpers.c)
    target_link_libraries(test_helpers flint)

    foreach (build_dir IN LISTS BUILD_DIRS CMAKE_SOURCE_DIR)
        file(GLOB TEST_FILES "${build_dir}/test/*.c")
        foreach(test_file IN LISTS TEST_FILES)
            file(RELATIVE_PATH test_name ${CMAKE_SOURCE_DIR} ${test_file})
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
                ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )
        endforeach()
    endforeach ()
endif()


if(BUILD_DOCS)
    find_package(Sphinx REQUIRED)
    file(GLOB DOC_SOURCES doc/source/*.rst)
    add_custom_target(html
        COMMAND ${SPHINX_EXECUTABLE} -b html "${CMAKE_SOURCE_DIR}/doc/source" "${CMAKE_BINARY_DIR}/html"
        SOURCES ${DOC_SOURCES})  
    add_custom_target(latex
        COMMAND ${SPHINX_EXECUTABLE} -b latex "${CMAKE_SOURCE_DIR}/doc/source" "${CMAKE_BINARY_DIR}/latex"
        SOURCES ${DOC_SOURCES})  
    add_custom_target(pdf DEPENDS latex COMMAND make WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/latex") 
endif()
