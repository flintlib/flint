################################################################################
# preamble
################################################################################

define(FLINT_COPYRIGHT,[[
Copyright (C) 2023 Albin Ahlbäck

This file is part of FLINT.

FLINT is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License (LGPL) as published
by the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.  See <https://www.gnu.org/licenses/>.

]])
define(FLINT_COPYRIGHT_C,[[/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/]])

AC_COPYRIGHT(FLINT_COPYRIGHT)
AH_TOP(FLINT_COPYRIGHT_C)

AC_PREREQ(2.59)

m4_define(FLINT_VERSION,m4_include([VERSION]))

AC_INIT(FLINT, FLINT_VERSION, [https://github.com/flintlib/flint/issues/], flint, [https://flintlib.org/])
AC_CONFIG_AUX_DIR([config])
AC_LANG(C)

# Every major version of Flint (1.0, 1.1, 1.2, etc.) gets
# a new major soname number.

# Every time a release happens where interfaces are added
# to Flint, the minor soname number gets incremented.

# Every patch that is made that doesn't change the Flint
# interface updates the patch number.

# However, if backwards incompatible changes are made, both
# the minor and patch numbers are set to 0. (Hopefully this
# doesn't happen in practice.)

# flint => soname
# 2.5.0 => 13.5.0
# 2.5.1 => 13.5.1
# 2.5.2 => 13.5.2
# 2.6.0 => 14.0.0
# 2.6.1 => 14.0.1
# 2.6.2 => 14.0.2
# 2.6.3 => 14.0.3
# 2.7.0 => 15.0.0
# 2.7.1 => 15.0.1
# 2.8.0 => 16.1.0
# 2.8.1 => 16.1.1
# 2.8.2 => 16.1.2
# 2.8.3 => 16.1.3
# 2.8.4 => 16.1.4
# 2.8.5 => 16.1.5
# 2.9.0 => 17.0.0
# 3.0.0 => 18.0.0

# NOTE: This must be after AC_INIT
FLINT_MAJOR_SO=18
FLINT_MINOR_SO=0
FLINT_PATCH_SO=0

# Make sure that we are in the right directory
AC_CONFIG_SRCDIR(src/fmpz/link/fmpz_reentrant.c)

################################################################################
# initialize libtool
################################################################################

# Later we want to overwrite the default CFLAGS. Check now whether they are set
# or not as Libtool overwrites CFLAGS.
if test -n "${CFLAGS+x}";
then
    cflags_set="yes"
else
    cflags_set="no"
fi

if test -n "${TESTCFLAGS+x}";
then
    testcflags_set="yes"
else
    testcflags_set="no"
fi

if test -n "${CXXFLAGS+x}";
then
    cxxflags_set="yes"
else
    cxxflags_set="no"
fi

LT_INIT

if test "$cflags_set" = "no";
then
    CFLAGS=
fi

if test "$cxxflags_set" = "no";
then
    CXXFLAGS=
fi

################################################################################
# build system
################################################################################

if test -n "$target_alias";
then
    AC_MSG_ERROR([FLINT does not support the use of specifying --target.])
fi

dnl Get system triplet
dnl NOTE: This is already invoked from LT_INIT
dnl AC_CANONICAL_HOST

################################################################################
# configure headers
################################################################################

AC_CONFIG_HEADERS(src/config.h src/flint-config.h)

################################################################################
# features
################################################################################

AC_ARG_ENABLE(pthread,
[AS_HELP_STRING([--enable-pthread],[Enable pthread [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-pthread. Need yes or no.])
    ;;
esac],
enable_pthread="yes")

AC_ARG_ENABLE(reentrant,
[AS_HELP_STRING([--enable-reentrant],[Build reentrant version of library [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-reentrant. Need yes or no.])
    ;;
esac],
enable_reentrant="no")

# Synonym for thread-safe. Only here for say to users that it is deprecated.
AC_ARG_ENABLE(tls)
AC_ARG_ENABLE(thread-safe,
[AS_HELP_STRING([--enable-thread-safe],[Enable thread-local storage [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-tls. Need yes or no.])
    ;;
esac],
enable_thread_safe="yes")

if test -n "$enable_tls";
then
    AC_MSG_ERROR([--enable-tls is an old feature. Please use --enable-thread-safe instead.])
fi

AC_ARG_ENABLE(assert,
[AS_HELP_STRING([--enable-assert],[Enable use of asserts [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-assert. Need yes or no.])
    ;;
esac],
enable_assert="no")

AC_ARG_ENABLE(coverage,
[AS_HELP_STRING([--enable-coverage],[Enable test coverage [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-coverage. Need yes or no.])
    ;;
esac],
enable_coverage="no")

AC_ARG_ENABLE(debug,
[AS_HELP_STRING([--enable-debug],[Compile FLINT with debug information [default=yes]])],
[case $enableval in
yes|no)
    if test "$enableval" = "no" && test "$enable_coverage" = "yes";
    then
        AC_MSG_ERROR([Debug information is required by coverage.])
    fi
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-debug. Need yes or no.])
    ;;
esac],
enable_debug="yes")

AC_ARG_ENABLE(dependency-tracking,
[AS_HELP_STRING([--enable-dependency-tracking],[Enable GCC automated dependency tracking [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-dependency-tracking. Need yes or no.])
    ;;
esac],
enable_dependency_tracking="yes")

AC_ARG_ENABLE(gmp-internals,
[AS_HELP_STRING([--enable-gmp-internals],[Enable calling GMP internals directly [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-gmp-internals. Need yes or no.])
    ;;
esac],
enable_gmp_internals="yes")

AC_ARG_ENABLE(avx2,
[AS_HELP_STRING([--enable-avx2],[Enable AVX2 instructions [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-avx2. Need yes or no.])
    ;;
esac],
enable_avx2="no")

AC_ARG_ENABLE(avx512,
[AS_HELP_STRING([--enable-avx512],[Enable AVX512 instructions [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-avx512. Need yes or no.])
    ;;
esac],
enable_avx512="no")

# Note: This is maintainer level only. Currently only used for Nemo CI.
AC_ARG_ENABLE(mpfr-check,[],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-mpfr-check. Need yes or no.])
    ;;
esac],
enable_mpfr_check="yes")

################################################################################
# packages
################################################################################

AC_ARG_WITH(gmp_include,
[AS_HELP_STRING([--with-gmp-include=DIR],[GMP include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
gmp_include_path="$withval"
)

AC_ARG_WITH(gmp_lib,
[AS_HELP_STRING([--with-gmp-lib=DIR],[GMP library directory])],
LDFLAGS="$LDFLAGS -L$withval"
gmp_lib_path="$withval"
)

AC_ARG_WITH(gmp,
[AS_HELP_STRING([--with-gmp[[=DIR]]],[GMP install directory])],
if test "$with_gmp" = "no";
then
    AC_MSG_FAILURE([The use of GMP is mandatory.])
elif test "$with_gmp" = "yes";
then
    :
elif test -z "$gmp_include_path" && test -z "$gmp_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    gmp_include_path="$withval/include"
    gmp_lib_path="$withval/lib"
    with_gmp="yes"
else
    AC_MSG_FAILURE([Cannot use --with-gmp along with --with-gmp-include or --with-gmp-lib.])
fi,
with_gmp="yes"
)

AC_ARG_WITH(mpfr_include,
[AS_HELP_STRING([--with-mpfr-include=DIR],[MPFR include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
mpfr_include_path="$withval"
)

AC_ARG_WITH(mpfr_lib,
[AS_HELP_STRING([--with-mpfr-lib=DIR],[MPFR library directory])],
LDFLAGS="$LDFLAGS -L$withval"
mpfr_lib_path="$withval"
)

AC_ARG_WITH(mpfr,
[AS_HELP_STRING([--with-mpfr[[=DIR]]],[MPFR install directory])],
if test "$with_mpfr" = "no";
then
    AC_MSG_FAILURE([The use of MPFR is mandatory.])
elif test "$with_mpfr" = "yes";
then
    :
elif test -z "$mpfr_include_path" && test -z "$mpfr_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    mpfr_include_path="$withval/include"
    mpfr_lib_path="$withval/lib"
    with_mpfr="yes"
else
    AC_MSG_FAILURE([Cannot use --with-mpfr along with --with-mpfr-include or --with-mpfr-lib.])
fi,
with_mpfr="yes"
)

AC_ARG_WITH(blas_include,
[AS_HELP_STRING([--with-blas-include=DIR],[Use BLAS and specify its include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
blas_include_path="$withval"
)

AC_ARG_WITH(blas_lib,
[AS_HELP_STRING([--with-blas-lib=DIR],[Use BLAS and specify its library directory])],
LDFLAGS="$LDFLAGS -L$withval"
blas_lib_path="$withval"
)

AC_ARG_WITH(blas,
[AS_HELP_STRING([--with-blas[[=DIR]]],[Use BLAS and specify its install directory])],
if test "$with_blas" = "yes";
then
    :
elif test "$with_blas" = "no";
then
    if test -n "$blas_include_path" || test -n "$blas_lib_path";
    then
        AC_MSG_FAILURE([--with-blas-include or --with-blas-lib was specified but --with-blas=no.])
    fi
elif test -z "$blas_include_path" && test -z "$blas_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    blas_include_path="$withval/include"
    blas_lib_path="$withval/lib"
    with_blas="yes"
else
    AC_MSG_FAILURE([Cannot use --with-blas along with --with-blas-include or --with-blas-lib.])
fi,
with_blas="no"
)

if test "$with_blas" != "no" || test -n "$blas_include_path" || test -n "$blas_lib_path";
then
    with_blas="yes"
fi

AC_ARG_WITH(gc_include,
[AS_HELP_STRING([--with-gc-include=DIR],[Use GC and specify its include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
gc_include_path="$withval",
with_gc_include="no"
)

AC_ARG_WITH(gc_lib,
[AS_HELP_STRING([--with-gc-lib=DIR],[Use GC and specify its library directory])],
LDFLAGS="$LDFLAGS -L$withval"
gc_lib_path="$withval",
with_gc_lib="no"
)

AC_ARG_WITH(gc,
[AS_HELP_STRING([--with-gc[[=DIR]]],[Use GC and optionally specify its install directory])],
if test "$with_gc" = "yes";
then
    :
elif test "$with_gc" = "no";
then
    if test -n "$gc_include_path" || test -n "$gc_lib_path";
    then
        AC_MSG_FAILURE([--with-gc-include or --with-gc-lib was specified but --with-gc=no.])
    fi
elif test -z "$gc_include_path" && test -z "$gc_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    gc_include_path="$withval/include"
    gc_lib_path="$withval/lib"
    with_gc="yes"
else
    AC_MSG_FAILURE([Cannot use --with-gc along with --with-gc-include or --with-gc-lib.])
fi,
with_gc="no"
)

if test "$with_gc" != "no" || test -n "$gc_include_path" || test -n "$gc_lib_path";
then
    with_gc="yes"
fi

AC_ARG_WITH(ntl_include,
[AS_HELP_STRING([--with-ntl-include=DIR],[Use NTL (in tests) and specify its include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
ntl_include_path="$withval"
)

AC_ARG_WITH(ntl_lib,
[AS_HELP_STRING([--with-ntl-lib=DIR],[Use NTL (in tests) and specify its library directory])],
LDFLAGS="$LDFLAGS -L$withval"
ntl_lib_path="$withval"
)

AC_ARG_WITH(ntl,
[AS_HELP_STRING([--with-ntl[[=DIR]]],[Use NTL (in tests) and optionally specify its install directory])],
if test "$with_ntl" = "yes";
then
    :
elif test "$with_ntl" = "no";
then
    if test -n "$ntl_include_path" || test -n "$ntl_lib_path";
    then
        AC_MSG_FAILURE([--with-ntl-include or --with-ntl-lib was specified but --with-ntl=no.])
    fi
elif test -z "$ntl_include_path" && test -z "$ntl_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    ntl_include_path="$withval/include"
    ntl_lib_path="$withval/lib"
    with_ntl="yes"
else
    AC_MSG_FAILURE([Cannot use --with-ntl along with --with-ntl-include or --with-ntl-lib.])
fi,
with_ntl="no"
)

if test "$with_ntl" != "no" || test -n "$ntl_include_path" || test -n "$ntl_lib_path";
then
    with_ntl="yes"
fi

################################################################################
# programs
################################################################################

dnl NOTE: The following line is invoked from LT_INIT
dnl AC_PROG_CC

if test "$with_ntl" = "yes";
then
    AC_PROG_CXX
    AC_LANG_PUSH(C++)
    if test "$cxxflags_set" = "no";
    then
        AX_CXX_CHECK_COMPILE_FLAG([-Wall],[CXXFLAGS="-Wall $CXXFLAGS"])
        AX_CXX_CHECK_COMPILE_FLAG([-Werror=implicit-function-declaration],[CXXFLAGS="-Werror=implicit-function-declaration $CXXFLAGS"])
        AX_CXX_CHECK_COMPILE_FLAG([-O2],[CXXFLAGS="-O2 $CXXFLAGS"])
        if test "$host" = "$build";
        then
            AX_CXX_CHECK_COMPILE_FLAG([-march=native],[CXXFLAGS="-march=native $CXXFLAGS"])
            fi
        AX_CXX_CHECK_COMPILE_FLAG([-std=c++11],[CXXFLAGS="-std=c++11 $CXXFLAGS"],
            AC_MSG_ERROR([Couldn't compile with the C++11 standard needed by the NTL interface.])
        )
        if test "$enable_debug" = "yes";
        then
            if "$ac_cv_prog_cxx_g" = "yes";
            then
                CXXFLAGS="-g $CXXFLAGS"
            else
                AC_MSG_ERROR([$CXX does not support the flag -g needed to generate debug information])
            fi
        fi
    fi
    AC_LANG_POP
fi

AC_PROG_MKDIR_P

################################################################################
# environment variables
################################################################################

AC_ARG_VAR([TESTCFLAGS], [Choose different C compiler flags for tests [default=CFLAGS]. Can be useful to compile library at a specific optimization level while compiling tests fast.])

AC_ARG_VAR(ABI, [Desired ABI])

AC_ARG_VAR(LDCONFIG, [ldconfig tool])

################################################################################
# check programs and system
################################################################################

AC_C_INLINE

if test "$compiler_c_o" = "no";
then
    AC_MSG_ERROR([Compiler must support both -c and -o simultaneously!])
fi

AC_C_BIGENDIAN(flint_big_endian=yes)

case "$host_os" in
    darwin*)
        flint_lib="libflint.dylib"
        flint_lib_major="libflint.$FLINT_MAJOR_SO.dylib"
        flint_lib_full="libflint.$FLINT_MAJOR_SO.$FLINT_MINOR_SO.dylib"
        extra_shared_flags="-install_name '`pwd`/$flint_lib_full' -compatibility_version $FLINT_MAJOR_SO.$FLINT_MINOR_SO -current_version $FLINT_MAJOR_SO.$FLINT_MINOR_SO.$FLINT_PATCH_SO"
        flint_dylib="yes"
        ;;
    cygwin|mingw*|msys)
        flint_lib="libflint.dll"
        flint_lib_major="libflint.$FLINT_MAJOR_SO.dll"
        flint_lib_full="libflint.$FLINT_MAJOR_SO.$FLINT_MINOR_SO.$FLINT_PATCH_SO.dll"
        flint_implib="$flint_lib.$libext"
        extra_shared_flags="-static-libgcc $wl--export-all-symbols,--out-implib,$flint_implib"
        flint_dlllib="yes"
        ;;
    *)
        flint_lib="libflint.so"
        flint_lib_major="libflint.so.$FLINT_MAJOR_SO"
        flint_lib_full="libflint.so.$FLINT_MAJOR_SO.$FLINT_MINOR_SO.$FLINT_PATCH_SO"
        extra_shared_flags="$wl-soname,$flint_lib_major"
        flint_solib="yes"
        ;;
esac

flint_lib_static="libflint.$libext"

if test "`uname -s`" = "android";
then
    extra_shared_flags="$wl-soname,$flint_lib"
fi

AC_SUBST(WL,$wl)

AC_SUBST(FLINT_LIB,$flint_lib)
AC_SUBST(FLINT_LIB_MAJOR,$flint_lib_major)
AC_SUBST(FLINT_LIB_FULL,$flint_lib_full)
AC_SUBST(FLINT_IMPLIB,$flint_implib)
AC_SUBST(FLINT_LIB_STATIC,$flint_lib_static)
AC_SUBST(EXTRA_SHARED_FLAGS,$extra_shared_flags)

if test "$flint_dlllib" = "yes";
then
    AC_SUBST(FLINT_DLLLIB,1)
else
    AC_SUBST(FLINT_DLLLIB,0)
fi

if test "$flint_dylib" = "yes";
then
    AC_SUBST(FLINT_DYLIB,1)
else
    AC_SUBST(FLINT_DYLIB,0)
fi

if test "$flint_solib" = "yes";
then
    AC_SUBST(FLINT_SOLIB,1)
else
    AC_SUBST(FLINT_SOLIB,0)
fi

case "$host_os" in
    cygwin|mingw*|msys)
        if test "$DLLTOOL" = "false";
        then
            AC_MSG_ERROR(["Couldn't find any dlltool that is required to create import libraries for
Windows-type systems."])
        fi
        ;;
    *)
        ;;
esac

case "$host_cpu" in
    x86_64|sparc64|powerpc64*)
        fft_tuning_in="src/fft_tuning64.in"
        ;;
    *)
        fft_tuning_in="src/fft_tuning32.in"
        ;;
esac
AC_CONFIG_FILES([src/fft_tuning.h:$fft_tuning_in],[],[fft_tuning_in="$fft_tuning_in"])
AC_SUBST(FFT_TUNING_IN, $fft_tuning_in)

if test -z "$LDCONFIG";
then
    AC_MSG_CHECKING(for ldconfig)
    case "$build_os" in
        mingw*|cygwin|msys|darwin*|freebsd*)
            LDCONFIG="true"
            ;;
        *)
            LDCONFIG="true"
            if test "$FLINT_SOLIB" = "1";
            then
                if test -x ldconfig > /dev/null 2>&1;
                then
                    LDCONFIG="ldconfig"
                elif test -x /sbin/ldconfig;
                then
                    LDCONFIG="/sbin/ldconfig"
                fi
            fi
            ;;
    esac
    if test "$LDCONFIG" = "true";
    then
        ldconfig_res="skipped"
    else
        ldconfig_res="$LDCONFIG"
    fi
    AC_MSG_RESULT($ldconfig_res)
fi
AC_SUBST(LDCONFIG)

dnl FIXME: Improve this check to include more processors.
AC_MSG_CHECKING([if memory is strongly-ordered])
case "$host_cpu" in
    x86_64|x86|i386|i586)
        flint_know_strong_order="yes"
        AC_MSG_RESULT([yes])
        ;;
    *)
        flint_know_strong_order="no"
        AC_MSG_RESULT([unsure])
        ;;
esac

if test "$flint_know_strong_order" = "yes";
then
    AC_DEFINE(FLINT_KNOW_STRONG_ORDER,1,[Define if system is strongly ordered])
fi

################################################################################
# check headers
################################################################################

# Mandatory headers
# The following headers are checked previously:
#   stdio.h stdlib.h stdint.h string.h
AC_CHECK_HEADERS([stdarg.h math.h float.h errno.h mpfr.h gmp.h],,
                 AC_MSG_ERROR([Could not find a mandatory header!]))

# Check that GMP >= 6.2.1
AC_MSG_CHECKING([if version of GMP is greater than 6.2.1])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
        [#include <gmp.h>],
        [#if (__GNU_MP_VERSION < 6) \
          || (__GNU_MP_VERSION == 6 && __GNU_MP_VERSION_MINOR < 2) \
          || (__GNU_MP_VERSION == 6 && __GNU_MP_VERSION_MINOR == 2 && __GNU_MP_VERSION_PATCHLEVEL < 1)
         # error GMP version 6.2.1 or later is required
         #endif]
    )],
    AC_MSG_RESULT([yes]),
    AC_MSG_RESULT([no])
    AC_MSG_ERROR(["FLINT requires GMP version 6.2.1 or later."]))

AC_MSG_CHECKING([if GMP defines mp_limb_t as unsigned long long int])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
        [#include <gmp.h>],
        [#if !defined(_LONG_LONG_LIMB)
         # error mp_limb_t != unsigned long long int
         #endif]
    )],
    AC_MSG_RESULT([yes])
    gmpcompat_h_in="gmpcompat-longlong.h.in",
    AC_MSG_RESULT([no])
    gmpcompat_h_in="gmpcompat.h.in")

AC_CONFIG_FILES([src/gmpcompat.h:src/$gmpcompat_h_in],[],[gmpcompat_h_in="$gmpcompat_h_in"])
AC_SUBST(GMPCOMPAT_H_IN, $gmpcompat_h_in)

# Check that MPFR >= 4.1.0
AC_MSG_CHECKING([if version of MPFR is greater than 4.1.0])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
        [#include <mpfr.h>],
        [#if (MPFR_VERSION_MAJOR < 4) \
          || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR < 1) \
          || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR == 1 && MPFR_VERSION_PATCHLEVEL < 0)
         # error MPFR version 4.1.0 or later is required
         #endif]
    )],
    AC_MSG_RESULT([yes]),
    AC_MSG_RESULT([no])
    AC_MSG_ERROR(["FLINT requires MPFR version 4.2.0 or later."]))

# Optional headers
# The following headers are checked previously:
#   unistd.h
AC_CHECK_HEADERS([fenv.h alloca.h malloc.h sys/param.h windows.h arm_neon.h x86intrin.h immintrin.h pthread_np.h])

if test "$enable_pthread" = "yes";
then
    AC_CHECK_HEADER([pthread.h],,
                    AC_MSG_ERROR([Could not find pthread.h!]))
fi

if test "$enable_assert" = "yes";
then
    AC_CHECK_HEADER([assert.h],,
                    AC_MSG_ERROR([Could not find assert.h!]))
fi

if test "$with_blas" = "yes";
then
    AC_CHECK_HEADER([cblas.h],,
                    AC_MSG_ERROR([Could not find cblas.h!]))
fi

if test "$with_gc" = "yes";
then
    AC_CHECK_HEADER([gc.h],,
                    AC_MSG_ERROR([Could not find gc.h!]))
fi

if test "$with_ntl" = "yes";
then
    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([NTL/ZZ.h],,
                    AC_MSG_ERROR([Could not find NTL/ZZ.h!]))
    AC_LANG_POP([C++])
fi

################################################################################
# check libraries
################################################################################

AC_SEARCH_LIBS([atan2],[m],
   [],
   [AC_MSG_ERROR(["The C math library was not found!"])])

AC_SEARCH_LIBS([__gmpz_init],[gmp],
   [],
   [AC_MSG_ERROR(["GMP library was not found.  If you indeed have GMP installed, please
submit a bug report to <https://github.com/flintlib/flint/issues/> so
that we can either fix the issue or give a more proper error message."])])

if test "$enable_gmp_internals" = "yes";
then
    AC_SEARCH_LIBS([__gmpn_gcd_11],[gmp],
        [],
        [AC_MSG_ERROR(["`mpn_gcd_11' was not found in the GMP library. It is needed to enable GMP
internals."])])

    AC_SEARCH_LIBS([__gmpn_div_q],[gmp],
        [],
        [AC_MSG_ERROR([`mpn_div_q' was not found in the GMP library. It is needed to enable GMP
internals.])])

    AC_SEARCH_LIBS([__gmpn_mul_basecase],[gmp],
        [],
        [AC_MSG_ERROR([`mpn_mul_basecase' was not found in the GMP library. It is needed to enable GMP
internals.])])

    AC_SEARCH_LIBS([__gmpn_modexact_1_odd],[gmp],
        [AC_DEFINE(FLINT_HAVE_MPN_MODEXACT_1_ODD,1,Define if system has mpn_modexact_1_odd)])
fi

if test "$enable_mpfr_check" = "yes";
then
    AC_SEARCH_LIBS([mpfr_init],[mpfr],
        [],
        [AC_MSG_ERROR([MPFR library was not found.  If you indeed have MPFR installed, please
submit a bug report to <https://github.com/flintlib/flint/issues/> so
that we can either fix the issue or give a more proper error message.])])

    AC_SEARCH_LIBS([mpfr_round_p],[mpfr],
        [],
        [AC_MSG_ERROR(["`mpfr_round_p' was not found in the MPFR library. It is required for FLINT."])])

    AC_SEARCH_LIBS([mpfr_mulhigh_n],[mpfr],
        [],
        [AC_MSG_ERROR(["`mpfr_mulhigh_n' was not found in the MPFR library. It is required for FLINT."])])

    AC_SEARCH_LIBS([mpfr_sqrhigh_n],[mpfr],
        [],
        [AC_MSG_ERROR(["`mpfr_sqrhigh_n' was not found in the MPFR library. It is required for FLINT."])])
else
    LIBS="$LIBS -lmpfr"
fi

if test "$enable_pthread" = "yes";
then
    # check if compiler accepts the preferred -pthread
    AX_CHECK_COMPILE_FLAG([-pthread],
        [LIBS="-pthread $LIBS"],
        # check if the compiler can do -lpthread
        AC_SEARCH_LIBS([pthread_create],[pthread],
            [],
            AC_MSG_ERROR([Couldn't link against POSIX threads via -pthread or -lpthread while
--enable-pthread was specified.])
        )
    )
fi

if test "$with_blas" = "yes";
then
    AC_SEARCH_LIBS([cblas_dgemm],[cblas openblas blas],[],
        [AC_MSG_ERROR([BLAS library was not found.  If you indeed have BLAS installed, please
submit a bug report to <https://github.com/flintlib/flint/issues/> so
that we can either fix the issue or give a more proper error message.])]
    )
fi

if test "$with_ntl" = "yes";
then
    AC_MSG_CHECKING([if linking simple code snippet of NTL works])
    SAVE_libs="$LIBS"
    LIBS="-lntl $LIBS"
    AC_LANG_PUSH([C++])
    dnl The following code snippet was taken from NTL's documentation
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([[#include <NTL/ZZ.h>]],
            [NTL::ZZ a, b, c;
std::cin >> a;
std::cin >> b;
c = (a+1)*(b+1);
std::cout << c << "\n";])],
        AC_MSG_RESULT([yes]),
        AC_MSG_RESULT([no])
        [AC_MSG_ERROR([NTL library was not found.  If you indeed have NTL installed, please
submit a bug report to <https://github.com/flintlib/flint/issues/> so
that we can either fix the issue or give a more proper error message.])]
    )
    AC_LANG_POP
    LIBS="$SAVE_libs"
fi

################################################################################
# check settings and environment
################################################################################

AC_SUBST(PIC_FLAG,$pic_flag)

AC_MSG_CHECKING([if cpu_set_t is supported])
AC_LINK_IFELSE([AC_LANG_PROGRAM([#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>],
[cpu_set_t s;CPU_ZERO(&s); pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), 0);])],
[AC_DEFINE(FLINT_USES_CPUSET,1,Define if system has cpu_set_t)
 AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])])

# ensure alloca works
AC_MSG_CHECKING([if alloca works])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#  ifdef _MSC_VER
#   include <malloc.h>
#   define alloca _alloca
#  else
#   include <alloca.h> /* We assume then that you have alloca.h */
#  endif
# endif
#endif],
[char * p = alloca(128);])],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])
AC_MSG_ERROR([Couldn't find alloca, which is required for FLINT. Please submit a bug
report to <https://github.com/flintlib/flint/issues/> and specify your
operating system.])])

################################################################################
# CFLAGS
################################################################################

save_CFLAGS="$CFLAGS"
CFLAGS=""

# The following is needed for Clang to check for unknown options.
AX_CHECK_COMPILE_FLAG([-Werror=unknown-warning-option], [CFLAGS="-Werror=unknown-warning-option"])

if test "$enable_coverage" = "yes";
then
    AX_CHECK_COMPILE_FLAG([--coverage],
        [save_CFLAGS="--coverage $save_CFLAGS"],
        [AC_MSG_ERROR([$CC does not to support test coverage flags])])
    LDFLAGS="--coverage $LDFLAGS"

    if test "$enable_pthread" = "yes";
    then
        AX_CHECK_COMPILE_FLAG([-fprofile-update=atomic],
            [save_CFLAGS="-fprofile-update=atomic $save_CFLAGS"],
            [AC_MSG_ERROR([$CC does not to support multi-threaded test coverage flags])])
    fi
fi

if test "$enable_debug" = "yes";
then
    if test "$ac_cv_prog_cc_g" = "yes";
    then
        save_CFLAGS="-g $save_CFLAGS"
    else
        AC_MSG_ERROR([$CC does not support the flag -g needed to generate debug information])
    fi
fi

if test -n "$ABI";
then
    AX_CHECK_COMPILE_FLAG([-m$ABI],
    [save_CFLAGS="-m$ABI $save_CFLAGS"],
    [AC_MSG_ERROR([$CC does not support the ABI flag -m$ABI])])
fi
AC_SUBST(ABI_FLAG)

if test "$enable_avx2" = "yes";
then
    AX_CHECK_COMPILE_FLAG([-mavx2],
    [save_CFLAGS="-mavx2 $save_CFLAGS"
     AC_DEFINE(FLINT_HAVE_AVX2,1,[Define if system has AVX2])],
    [AC_MSG_ERROR([$CC does not support the flag -mavx2 needed for AVX2 instructions])])
    AX_CHECK_COMPILE_FLAG([-mfma], [save_CFLAGS="-mfma $save_CFLAGS"])
fi

if test "$enable_avx512" = "yes";
then
    AX_CHECK_COMPILE_FLAG([-mavx512f],
    [save_CFLAGS="-mavx512f $save_CFLAGS"
     AC_DEFINE(FLINT_HAVE_AVX512,1,[Define if system has AVX512])],
    [AC_MSG_ERROR([$CC does not support the flag -mavx512f needed for AVX512 instructions])])
fi

case $host_cpu in
    sparc*)
        AX_CHECK_COMPILE_FLAG(-mno-relax,[DEFAULT_CFLAGS="-mno-relax $DEFAULT_CFLAGS"])
        ;;
    *)
        ;;
esac

AX_CHECK_COMPILE_FLAG([-Wall],[DEFAULT_CFLAGS="-Wall $DEFAULT_CFLAGS"])
AX_CHECK_COMPILE_FLAG([-Wno-stringop-overread],[DEFAULT_CFLAGS="-Wno-stringop-overread $DEFAULT_CFLAGS"])
AX_CHECK_COMPILE_FLAG([-Wno-stringop-overflow],[DEFAULT_CFLAGS="-Wno-stringop-overflow $DEFAULT_CFLAGS"])
AX_CHECK_COMPILE_FLAG([-Werror=implicit-function-declaration],[DEFAULT_CFLAGS="-Werror=implicit-function-declaration $DEFAULT_CFLAGS"])
AX_CHECK_COMPILE_FLAG([-O3],[DEFAULT_CFLAGS="-O3 $DEFAULT_CFLAGS"])
if test "$host" = "$build";
then
    AX_CHECK_COMPILE_FLAG([-march=native],[DEFAULT_CFLAGS="-march=native $DEFAULT_CFLAGS"])
fi
AX_CHECK_COMPILE_FLAG([-std=c11],[DEFAULT_CFLAGS="-std=c11 $DEFAULT_CFLAGS"])
AX_CHECK_COMPILE_FLAG([-pedantic],[DEFAULT_CFLAGS="-pedantic $DEFAULT_CFLAGS"])

dnl NOTE: The following flags shouldn't be default. They are leftovers from
dnl previous configure script
AX_CHECK_COMPILE_FLAG([-funroll-loops],[DEFAULT_CFLAGS="-funroll-loops $DEFAULT_CFLAGS"])

if test "$cflags_set" = "no";
then
    CFLAGS="$save_CFLAGS $DEFAULT_CFLAGS"
else
    CFLAGS="$save_CFLAGS"
    AX_CHECK_COMPILE_FLAG([$CFLAGS],[],AC_MSG_ERROR(["Couldn't compile with given CFLAGS!"]))
fi

if test "$testcflags_set" = "no";
then
    TESTCFLAGS="$CFLAGS"
else
    TESTCFLAGS="$TESTCFLAGS $save_CFLAGS"
fi

CPPFLAGS="-I./src $CPPFLAGS"

################################################################################
# fft_small module
################################################################################

enable_fft_small="no"

AC_MSG_CHECKING([whether fft_small module is available])

if test "$ac_cv_header_arm_neon_h" = "yes";
then
    AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([],[
#ifndef __ARM_NEON
#error Dead man
error
#endif
         ])],
        [AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([],[
#ifndef __GNUC__
#error Dead man
error
#endif
             ])],
            [enable_fft_small="yes"]
        )]
    )
fi

if test "$ac_cv_header_x86intrin_h" = "yes" && test "$ac_cv_header_immintrin_h" = "yes";
then
    AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([],[
#ifndef __AVX2__
#error Dead man
error
#endif
         ])],
        [AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([],[
#ifndef __GNUC__
#error Dead man
error
#endif
             ])],
            [AC_COMPILE_IFELSE(
                [AC_LANG_PROGRAM([#include <gmp.h>],[
#if GMP_LIMB_BITS == 32
#error Dead man
error
#endif
                 ])],
                [enable_fft_small="yes"]
            )]
        )]
    )
fi

AC_MSG_RESULT([$enable_fft_small])

if test "$enable_fft_small" = "yes";
then
    AC_SUBST(FFT_SMALL, [fft_small\ \ \ ])
    AC_DEFINE(FLINT_HAVE_FFT_SMALL, 1, [Define to use the fft_small module])
else
    AC_SUBST(FFT_SMALL, [\ \ \ \ \ \ \ \ \ \ \ \ ])
fi

################################################################################
# substitutions and definitions
################################################################################

AC_SUBST(FLINT_VERSION_FULL, FLINT_VERSION)

# split version into major/minor/patch using POSIX variable substitutions
tail=FLINT_VERSION
FLINT_MAJOR=${tail%%.*}
tail=${tail#*.}
FLINT_MINOR=${tail%%.*}
tail=${tail#*.}
FLINT_PATCH=${tail%-*}

AC_SUBST(FLINT_MAJOR)
AC_SUBST(FLINT_MINOR)
AC_SUBST(FLINT_PATCH)
AC_SUBST(FLINT_MAJOR_SO)
AC_SUBST(FLINT_MINOR_SO)
AC_SUBST(FLINT_PATCH_SO)

AC_SUBST(GMP_LIB_PATH,$gmp_lib_path)
AC_SUBST(MPFR_LIB_PATH,$mpfr_lib_path)
AC_SUBST(BLAS_LIB_PATH,$blas_lib_path)
AC_SUBST(GC_LIB_PATH,$gc_lib_path)
AC_SUBST(NTL_LIB_PATH,$ntl_lib_path)

AC_SUBST(LIB_CPPFLAGS,["-DFLINT_NOSTDIO -DFLINT_NOSTDARG"])

# Configure fmpz's memory manager
if test "$with_gc" = "yes";
then
    if test "$enable_thread_safe" = "yes";
    then
        AC_MSG_ERROR(["The Boehm-Demers-Weise garbage collector does not support thread-local storage!"])
    fi
    fmpz_c="src/fmpz/link/fmpz_gc.c"
else
    if test "$enable_reentrant" = "yes";
    then
        fmpz_c="src/fmpz/link/fmpz_reentrant.c"
    else
        fmpz_c="src/fmpz/link/fmpz_single.c"
    fi
fi

AC_CONFIG_FILES([src/fmpz/fmpz.c:$fmpz_c],[],[fmpz_c="$fmpz_c"])
AC_SUBST(FMPZ_C_IN, $fmpz_c)

if test "$enable_shared" = "yes";
then
    AC_SUBST(SHARED,1)
else
    AC_SUBST(SHARED,0)
fi

if test "$enable_static" = "yes";
then
    AC_SUBST(STATIC,1)
else
    AC_SUBST(STATIC,0)
fi

if test "$enable_coverage" = "yes";
then
    AC_SUBST(COVERAGE,1)
    AC_DEFINE(FLINT_COVERAGE,1,[Define to enable coverage.])
else
    AC_SUBST(COVERAGE,0)
fi

if test "$enable_dependency_tracking" = "yes";
then
    AC_SUBST(WANT_DEPS,1)
else
    AC_SUBST(WANT_DEPS,0)
fi

if test "$with_ntl" = "yes";
then
    AC_SUBST(WANT_NTL,1)
else
    AC_SUBST(WANT_NTL,0)
fi

if test "$with_blas" = "yes";
then
    AC_DEFINE(FLINT_USES_BLAS,1,[Define to enable BLAS.])
fi

if test "$with_gc" = "yes";
then
    AC_DEFINE(FLINT_USES_GC,1,[Define to enable the Boehm-Demers-Weise garbage collector.])
fi

if test "$flint_big_endian" = "yes";
then
    AC_DEFINE(FLINT_BIG_ENDIAN,1,[Define if system is big endian.])
fi

if test "$enable_pthread" = "yes";
then
    AC_DEFINE(FLINT_USES_PTHREAD,1,[Define to enable the use of pthread.])
fi

if test "$enable_thread_safe" = "yes";
then
    AC_DEFINE(FLINT_USES_TLS,1,[Define to enable thread-local storage.])
fi

if test "$enable_reentrant" = "yes";
then
    AC_DEFINE(FLINT_REENTRANT,1,[Define to enable reentrant.])
fi

if test "$enable_assert" = "yes";
then
    AC_DEFINE(FLINT_WANT_ASSERT,1,[Define to enable use of asserts.])
fi

if test "$enable_gmp_internals" = "yes";
then
    AC_DEFINE(FLINT_WANT_GMP_INTERNALS,1,[Define to enable use of GMP internals.])
fi

################################################################################
# epilog
################################################################################

AC_CONFIG_FILES([Makefile flint.pc src/flint.h])
AC_OUTPUT

dnl Shorten the original help message
m4_cleardivert([HELP_BEGIN])dnl
m4_divert_push([HELP_BEGIN])dnl

if test "$ac_init_help" = "long";
then
    cat <<_ACEOF
\`configure' configures m4_ifset([AC_PACKAGE_STRING],
			[AC_PACKAGE_STRING],
			[this package]) to adapt to many kinds of systems.

Usage: $[0] [[OPTION]]... [[VAR=VALUE]]...

[To assign environment variables (e.g., CC, CFLAGS...), specify them as
VAR=VALUE.  See below for descriptions of some of the useful variables.

Defaults for the options are specified in brackets.

Configuration:
  -h, --help              display this help and exit
      --help=short        display options specific to this package
      --help=recursive    display the short help of all the included packages
  -V, --version           display version information and exit
  -q, --quiet, --silent   do not print \`checking ...' messages
      --cache-file=FILE   cache test results in FILE [disabled]
  -C, --config-cache      alias for \`--cache-file=config.cache'
  -n, --no-create         do not create output files
      --srcdir=DIR        find the sources in DIR [configure dir or \`..']

Installation directories:
]AS_HELP_STRING([--prefix=PREFIX],
  [install architecture-independent files in PREFIX [$ac_default_prefix]])
AS_HELP_STRING([--exec-prefix=EPREFIX],
  [install architecture-dependent files in EPREFIX [PREFIX]])[

By default, \`make install' will install all the files in
\`$ac_default_prefix/include' and \`$ac_default_prefix/lib' etc.
You can specify an installation prefix other than \`$ac_default_prefix'
using \`--prefix', for instance \`--prefix=\$HOME'.

For better control, use the options below.

Fine tuning of the installation directories:
  --bindir=DIR            user executables [EPREFIX/bin]
  --libdir=DIR            object code libraries [EPREFIX/lib]
  --includedir=DIR        C header files [PREFIX/include]
_ACEOF

    cat <<\_ACEOF]
m4_divert_pop([HELP_BEGIN])dnl

dnl NOTE: There is supposed to be a missing `fi' here.
