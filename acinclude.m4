# ===========================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_check_compile_flag.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_COMPILE_FLAG(FLAG, [ACTION-SUCCESS], [ACTION-FAILURE], [EXTRA-FLAGS], [INPUT])
#
# DESCRIPTION
#
#   Check whether the given FLAG works with the current language's compiler
#   or gives an error.  (Warnings, however, are ignored)
#
#   ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on
#   success/failure.
#
#   If EXTRA-FLAGS is defined, it is added to the current language's default
#   flags (e.g. CFLAGS) when the check is done.  The check is thus made with
#   the flags: "CFLAGS EXTRA-FLAGS FLAG".  This can for example be used to
#   force the compiler to issue an error when a bad flag is given.
#
#   INPUT gives an alternative input source to AC_COMPILE_IFELSE.
#
#   NOTE: Implementation based on AX_CFLAGS_GCC_OPTION. Please keep this
#   macro in sync with AX_CHECK_{PREPROC,LINK}_FLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#   Copyright (c) 2023 Albin Ahlbäck
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 6

# This code is modified to print CC instead of *language* compiler

AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether $CC accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([m4_default([$5],[AC_LANG_PROGRAM()])],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_VAR_IF(CACHEVAR,yes,
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

# Test for CXX compiler
AC_DEFUN([AX_CXX_CHECK_COMPILE_FLAG],
[AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether $CXX accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([m4_default([$5],[AC_LANG_PROGRAM()])],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_VAR_IF(CACHEVAR,yes,
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CXX_CHECK_COMPILE_FLAGS


# Copyright (C) 1996-2024 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.


dnl  AX_INIT
dnl  -----------------------
dnl  If build directory is not source directory, this function throws if source
dnl  directory is already configured.

AC_DEFUN([AX_INIT],[dnl
if test "$ac_abs_confdir" != "`pwd`"; dnl ' Vim syntax fix
then
  if test -f $srcdir/config.status;
  then
    AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
  fi
fi])


dnl Copyright (C) 2024 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.

dnl NOTE: The first two patterns are taken from GMP. The license is down below.

define(X86_PATTERN,
[[i?86*-*-* | k[5-8]*-*-* | pentium*-*-* | athlon-*-* | viac3*-*-* | geode*-*-* | atom-*-*]])

define(X86_64_PATTERN,
[[athlon64-*-* | k8-*-* | k10-*-* | bobcat-*-* | jaguar*-*-* | bulldozer*-*-* | piledriver*-*-* | steamroller*-*-* | excavator*-*-* | zen*-*-* | pentium4-*-* | atom-*-* | silvermont-*-* | goldmont-*-* | tremont-*-* | core2-*-* | corei*-*-* | x86_64-*-* | nano-*-* | nehalem*-*-* | westmere*-*-* | sandybridge*-*-* | ivybridge*-*-* | haswell*-*-* | broadwell*-*-* | skylake*-*-* | kabylake*-*-* | icelake*-*-* | tigerlake*-*-* | rocketlake*-*-* | alderlake*-*-* | raptorlake*-*-* | x86_64v[1234]-*-*]])

define(X86_64_ADX_PATTERN,
[[zen[1234]-*-* | coreibwl-*-* | broadwell-*-* | skylake-*-* | skylake_server-*-* | cannonlake-*-* | kabylake-*-* | icelake-*-* | icelake_server-*-* | rocketlake-*-* | tigerlake-*-* | alderlake-*-* | raptorlake-*-* | knightslanding-*-* | sapphirerapids-*-* | cometlake-*-*]])

define(ARM64_PATTERN,
[[armcortexa53-*-* | armcortexa53neon-*-* | armcortexa55-*-* | armcortexa55neon-*-* | armcortexa57-*-* | armcortexa57neon-*-* | armcortexa7[2-9]-*-* | armcortexa7[2-9]neon-*-* | armexynosm1-*-* | armthunderx-*-* | armxgene1-*-* | aarch64*-*-* | applem[1-9]*-*-* | armv8*-*-*]])

define(SLOW_VROUNDPD_PATTERN,
[[haswell* | broadwell* | skylake* | kabylake* | icelake* | tigerlake* | rocketlake* | alderlake* | raptorlake*]])

define(FAST_VROUNDPD_PATTERN,
[[znver[2-4]* | sandybridge* | ivybridge*]])



dnl  FLINT_SET_LIBFLAGS(lib,lib-path,include-path,[library_alias])
dnl  -----------------------
dnl  Sets lib_LDFLAGS and lib_CPPFLAGS to include and link library.
dnl  If pkg-config is available, lib_CFLAGS, lib_LIBS, lib_libdir,
dnl  lib_includedir are also set and also sets lib_LIBS with the appropriate
dnl  `-l' flag(s).  Else, straight up use lib-path and include-path.

AC_DEFUN([FLINT_SET_LIBFLAGS],
[tmpalias=m4_default([$4],[$1])

if test "x$2" != "x";
then
    pkgconfigpath="PKG_CONFIG_PATH=$2/pkgconfig"
fi

if test "$enable_pkg_config" = "yes" && eval "$pkgconfigpath $PKG_CONFIG --path $1"; test "$?" = "0";
then
    # libdir
    tmp=`eval "$pkgconfigpath $PKG_CONFIG --variable=libdir $1"` dnl ' Fix Vim syntax
    eval ${tmpalias}_libdir="\${tmp}"
    echo "libdir = $tmp"

    # includedir
    tmp=`eval "$pkgconfigpath $PKG_CONFIG --variable=includedir $1"` dnl ' Fix Vim syntax
    eval ${tmpalias}_includedir="\${tmp}"
    echo "includedir = $tmp"

    # LIBS
    tmp=`eval "$pkgconfigpath $PKG_CONFIG --libs-only-l $1"` dnl ' Fix Vim syntax
    eval ${tmpalias}_LIBS="\${tmp}"
    echo "LIBS = $tmp"

    # LDFLAGS
    tmp=`eval "$pkgconfigpath $PKG_CONFIG --libs-only-L $1"` dnl ' Fix Vim syntax
    eval ${tmpalias}_LDFLAGS="\${tmp}"
    echo "LDFLAGS = $tmp"

    # CPPFLAGS
    tmp=`eval "$pkgconfigpath $PKG_CONFIG --cflags-only-other $1" | sed -n 's/\(-D\w\+\)\(\|=\w\+\)/\n\1\2\n/gp' | sed -n '/^-D/p'` dnl ' Fix Vim syntax
    tmp2=`eval "$pkgconfigpath $PKG_CONFIG --cflags-only-I $1"` dnl ' Fix Vim syntax
    tmp="$tmp $tmp2"
    eval ${tmpalias}_CPPFLAGS="\${tmp}"
    echo "CPPFLAGS = $tmp"

    # CFLAGS
    tmp=`eval "$pkgconfigpath $PKG_CONFIG --cflags-only-other $withpath $1" | sed 's/\(-D\w\+\)\(\|=\w\+\)//g' | sed 's/  / /g'` dnl ' Fix Vim syntax
    eval ${tmpalias}_CFLAGS="\${tmp}"
    echo "CFLAGS = $tmp"
else
    if test "x${2}" != "x";
    then
        eval ${tmpalias}_LDFLAGS="-L\${2}"
    fi

    if test "x${3}" != "x";
    then
        eval ${tmpalias}_CPPFLAGS="-I\${3}"
    fi
fi
])



dnl  FLINT_CC_IS_GCC([action-if-true],[action-if-false])
dnl  -----------------------
dnl  Checks if compiler is GCC.

AC_DEFUN([FLINT_CC_IS_GCC],
[AC_CACHE_CHECK([if compiler is GCC],
                flint_cv_cc_is_gcc,
[flint_cv_cc_is_gcc="no"
AC_PREPROC_IFELSE([AC_LANG_PROGRAM([
#if !(defined(__GNUC__) && !defined(__clang__))
#error
error
#endif
],[])],
[flint_cv_cc_is_gcc="yes"])
])
AS_VAR_IF([flint_cv_cc_is_gcc],"yes",
    [m4_default([$1], :)],
    [m4_default([$2], :)])
])


dnl  FLINT_CC_IS_CLANG([action-if-true],[action-if-false])
dnl  -----------------------
dnl  Checks if compiler is clang.

AC_DEFUN([FLINT_CC_IS_CLANG],
[AC_CACHE_CHECK([if compiler is Clang],
                flint_cv_cc_is_clang,
[flint_cv_cc_is_clang="no"
AC_PREPROC_IFELSE([AC_LANG_PROGRAM([
#ifndef __clang__
#error
error
#endif
],[])],
[flint_cv_cc_is_clang="yes"])
])
AS_VAR_IF([flint_cv_cc_is_clang],"yes",
    [m4_default([$1], :)],
    [m4_default([$2], :)])
])


dnl  FLINT_CHECK_CPU_SET_T([action-if-true],[action-if-false])
dnl  -----------------------
dnl  Checks if cpu_set_t is supported.
dnl
dnl  FIXME: Does this cover all BSD systems?

AC_DEFUN([FLINT_CHECK_CPU_SET_T],
[AC_CACHE_CHECK([if cpu_set_t is supported],
                flint_cv_check_cpu_set_t,
[flint_cv_check_cpu_set_t="no"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>
],
[cpu_set_t s;
CPU_ZERO(&s);
pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), 0);])],
[flint_cv_check_cpu_set_t="yes"])
])
AS_VAR_IF([flint_cv_check_cpu_set_t],"yes",
    [m4_default([$1], :)],
    [m4_default([$2], :)])
])


dnl  FLINT_CHECK_NTL([action-if-true],[action-if-false])
dnl  -----------------------
dnl  Checks if linking with NTL works.

AC_DEFUN([FLINT_CHECK_NTL],
[AC_REQUIRE([AC_PROG_CXX])
AC_CACHE_CHECK([if linking with NTL works],
                flint_cv_check_ntl,
[flint_cv_check_ntl="no"
save_LIBS="$LIBS"
LIBS="-lntl $LIBS"
AC_LANG_PUSH([C++])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[#include <NTL/ZZ.h>
    ]], [NTL::ZZ a, b, c;
    std::cin >> a;
    std::cin >> b;
    c = (a+1)*(b+1);
    std::cout << c << "\n";])],
    [flint_cv_check_ntl="yes"])
AC_LANG_POP([C++])
LIBS="$save_LIBS"
])
AS_VAR_IF([flint_cv_check_ntl],"yes",
    [m4_default([$1], :)],
    [m4_default([$2], :)])
])


dnl  FLINT_PREPROC_IFELSE(input,[action-if-true],[action-if-false])
dnl  -----------------------
dnl  Runs preprocessor with CFLAGS.
dnl
dnl  FIXME: Autoconf states that some compilers do not accept CFLAGS in the
dnl  preprocessor. Which ones is it? GCC and Clang certainly does. If it does
dnl  not allow CFLAGS into preprocessor, we should just compile normally.

AC_DEFUN([FLINT_PREPROC_IFELSE],
[AC_REQUIRE([AC_PROG_CPP])
save_ac_cpp="$ac_cpp"
ac_cpp="$CPP $CFLAGS $CPPFLAGS"
AC_PREPROC_IFELSE($@)
ac_cpp="$ac_cpp"
])


dnl  FLINT_CHECK_GMP_H(MAJOR, MINOR, PATCHLEVEL)
dnl  -----------------------
dnl  Checks that gmp.h can be found and that its version fullfills the version
dnl  requirement.

AC_DEFUN([FLINT_CHECK_GMP_H],
[AC_CHECK_HEADER([gmp.h],,AC_MSG_ERROR([Could not find gmp.h]))
AC_MSG_CHECKING([if version of GMP is greater than $1.$2.$3])
AC_PREPROC_IFELSE([AC_LANG_PROGRAM(
        [#include <gmp.h>
        ],[#if (__GNU_MP_VERSION < $1) \
          || (__GNU_MP_VERSION == $1 && __GNU_MP_VERSION_MINOR < $2) \
          || (__GNU_MP_VERSION == $1 && __GNU_MP_VERSION_MINOR == $2 && __GNU_MP_VERSION_PATCHLEVEL < $3)
        # error GMP version $1.$2.$3 or later is required
        #endif
        ]
    )],
    AC_MSG_RESULT([yes]),
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([GMP version $1.$2.$3 or later is required.])
)
])


dnl  FLINT_CHECK_MPFR_H(MAJOR, MINOR, PATCHLEVEL)
dnl  -----------------------
dnl  Checks that mpfr.h can be found and that its version fullfills the version
dnl  requirement.

AC_DEFUN([FLINT_CHECK_MPFR_H],
[AC_REQUIRE([FLINT_CHECK_GMP_H])
AC_CHECK_HEADER([mpfr.h],,AC_MSG_ERROR([Could not find mpfr.h]))
AC_MSG_CHECKING([if version of MPFR is greater than $1.$2.$3])
AC_PREPROC_IFELSE([AC_LANG_PROGRAM(
        [#include <mpfr.h>
        ],[#if (MPFR_VERSION_MAJOR < $1) \
         || (MPFR_VERSION_MAJOR == $1 && MPFR_VERSION_MINOR < $2) \
         || (MPFR_VERSION_MAJOR == $1 && MPFR_VERSION_MINOR == $2 && MPFR_VERSION_PATCHLEVEL < $3)
        # error MPFR version $1.$2.$3 or later is required
        #endif
        ]
    )],
    AC_MSG_RESULT([yes]),
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([MPFR version $1.$2.$3 or later is required.])
)
])


dnl  FLINT_GMP_LONG_LONG_LIMB([action-success][,action-fail])
dnl  -----------------------
dnl  Check if GMP uses long long limb.

AC_DEFUN([FLINT_GMP_LONG_LONG_LIMB],
[AC_REQUIRE([FLINT_CHECK_GMP_H])
AC_CACHE_CHECK([if GMP defines mp_limb_t as unsigned long long int],
                flint_cv_gmp_long_long_limb,
[AC_PREPROC_IFELSE([AC_LANG_PROGRAM(
        [
        #include <gmp.h>
        ],[
        #if !defined(_LONG_LONG_LIMB)
        # error mp_limb_t != unsigned long long int
        #endif
        ]
    )],
    flint_cv_gmp_long_long_limb="yes",
    flint_cv_gmp_long_long_limb="no"
)])

AS_VAR_IF([flint_cv_gmp_long_long_limb],"yes",
    [m4_default([$1], :)],
    [m4_default([$2], :)])
])


dnl  FLINT_ABI
dnl  -----------------------
dnl  Checks what ABI to use. First checks if ABI is given. If none where given,
dnl  check GMP's configuration.

AC_DEFUN([FLINT_ABI],
[AC_REQUIRE([FLINT_CHECK_GMP_H])
AC_ARG_VAR(ABI, [Desired ABI])
AC_CACHE_CHECK([for desired ABI],
                flint_cv_abi,
[if test -n "$ABI";
then
    flint_cv_abi="$ABI"
else
    AC_PREPROC_IFELSE([AC_LANG_PROGRAM(
            [
            #include <gmp.h>
            ],[
            #if GMP_LIMB_BITS == 32
            #error Dead man
            error
            #endif
            ]
        )],
        flint_cv_abi="64",
        flint_cv_abi="32"
    )
fi
])
])


dnl  FLINT_HAVE_FFT_SMALL_ARM_H
dnl  -----------------------
dnl  Checks if system have headers for fft_small on Arm. Will only run if on
dnl  arm64.

AC_DEFUN([FLINT_HAVE_FFT_SMALL_ARM_H],
[case $host in
    ARM64_PATTERN)
        AC_CHECK_HEADERS([arm_neon.h],
            flint_cv_have_fft_small_arm_h="yes",
            flint_cv_have_fft_small_arm_h="no")
        ;;
esac])


dnl  FLINT_HAVE_FFT_SMALL_ARM_I
dnl  -----------------------
dnl  Checks if system supports Arm NEON instructions.

AC_DEFUN([FLINT_HAVE_FFT_SMALL_ARM_I],
[case $host in
    ARM64_PATTERN)
        AC_CACHE_CHECK([if system have required ARM instruction set for fft_small],
                        flint_cv_have_fft_small_arm_i,
            [FLINT_PREPROC_IFELSE([AC_LANG_SOURCE([
                    #if !defined(__ARM_NEON)
                    #error Dead man
                    error
                    #endif
                ])],
                flint_cv_have_fft_small_arm_i="yes",
                flint_cv_have_fft_small_arm_i="no")])
        ;;
esac])


dnl  FLINT_HAVE_FFT_SMALL_X86_H
dnl  -----------------------
dnl  Checks if system have headers for fft_small on x86.

AC_DEFUN([FLINT_HAVE_FFT_SMALL_X86_H],
[case $host in
    X86_64_PATTERN)
        AC_CHECK_HEADERS([immintrin.h],
            flint_cv_have_fft_small_x86_h="yes",
            flint_cv_have_fft_small_x86_h="no")
        ;;
esac])


dnl  FLINT_HAVE_FFT_SMALL_X86_I
dnl  -----------------------
dnl  Checks if system supports AVX2 instructions.

AC_DEFUN([FLINT_HAVE_FFT_SMALL_X86_I],
[case $host in
    X86_64_PATTERN)
        AC_CACHE_CHECK([if system have required x86_64 instruction set for fft_small],
                        flint_cv_have_fft_small_x86_i,
            [FLINT_PREPROC_IFELSE([AC_LANG_SOURCE([
                    #if !defined(__AVX2__)
                    #error Dead man
                    error
                    #endif
                ])],
                flint_cv_have_fft_small_x86_i="yes",
                flint_cv_have_fft_small_x86_i="no")])
        ;;
esac])


dnl  FLINT_CHECK_FFT_SMALL([action-success][,action-fail])
dnl  -----------------------
dnl  Checks if fft_small module is available.
dnl  Do "action-success" if this succeeds, "action-fail" if not.

AC_DEFUN([FLINT_CHECK_FFT_SMALL],
[AC_REQUIRE([FLINT_ABI])
AC_REQUIRE([FLINT_HAVE_FFT_SMALL_ARM_H])
AC_REQUIRE([FLINT_HAVE_FFT_SMALL_ARM_I])
AC_REQUIRE([FLINT_HAVE_FFT_SMALL_X86_H])
AC_REQUIRE([FLINT_HAVE_FFT_SMALL_X86_I])

AC_CACHE_CHECK([if system can use FLINT's fft_small module],
                flint_cv_check_fft_small,
[flint_cv_check_fft_small="no"
if test "$flint_cv_abi" = "64";
then
    if test "$flint_cv_have_fft_small_arm_h" = "yes" && test "$flint_cv_have_fft_small_arm_i" = "yes";
    then
        flint_cv_check_fft_small="yes"
    fi
    if test "$flint_cv_have_fft_small_x86_h" = "yes" && test "$flint_cv_have_fft_small_x86_i" = "yes";
    then
        flint_cv_check_fft_small="yes"
    fi
fi
])

AS_VAR_IF([flint_cv_check_fft_small],"yes",
    [m4_default([$1], :)],
    [m4_default([$2], :)])
])


dnl  GMP specific autoconf macros
dnl  (Taken from GMP 6.3.0)
dnl
dnl  Copyright 2000-2006, 2009, 2011, 2013-2018 Free Software Foundation, Inc.
dnl
dnl  This file is part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can redistribute it and/or modify
dnl  it under the terms of either:
dnl
dnl    * the GNU Lesser General Public License as published by the Free
dnl      Software Foundation; either version 3 of the License, or (at your
dnl      option) any later version.
dnl
dnl  or
dnl
dnl    * the GNU General Public License as published by the Free Software
dnl      Foundation; either version 2 of the License, or (at your option) any
dnl      later version.
dnl
dnl  or both in parallel, as here.
dnl
dnl  The GNU MP Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl  for more details.
dnl
dnl  You should have received copies of the GNU General Public License and the
dnl  GNU Lesser General Public License along with the GNU MP Library.  If not,
dnl  see https://www.gnu.org/licenses/.

dnl  GMP_INIT([M4-DEF-FILE])
dnl  -----------------------
dnl  Initializations for GMP config.m4 generation.
dnl
dnl  FIXME: The generated config.m4 doesn't get recreated by config.status.
dnl  Maybe the relevant "echo"s should go through AC_CONFIG_COMMANDS.

dnl NOTE: This is different from GMP.

AC_DEFUN([GMP_INIT],
[ifelse([$1], , gmp_configm4=config.m4, gmp_configm4="[$1]")
gmp_tmpconfigm4=cnfm4.tmp
gmp_tmpconfigm4i=cnfm4i.tmp
gmp_tmpconfigm4p=cnfm4p.tmp
rm -f $gmp_tmpconfigm4 $gmp_tmpconfigm4i $gmp_tmpconfigm4p

echo ["define(<CONFIG_TOP_SRCDIR>,<\`$srcdir'>)"] >>$gmp_tmpconfigm4

echo ["include][(CONFIG_TOP_SRCDIR\`/src/mpn_extras/asm-defs.m4')"] >>$gmp_tmpconfigm4i
])


dnl  GMP_FINISH
dnl  ----------
dnl  Create config.m4 from its accumulated parts.
dnl
dnl  __CONFIG_M4_INCLUDED__ is used so that a second or subsequent include
dnl  of config.m4 is harmless.
dnl
dnl  A separate ifdef on the angle bracket quoted part ensures the quoting
dnl  style there is respected.  The basic defines from gmp_tmpconfigm4 are
dnl  fully quoted but are still put under an ifdef in case any have been
dnl  redefined by one of the m4 include files.
dnl
dnl  Doing a big ifdef within asm-defs.m4 and/or other macro files wouldn't
dnl  work, since it'd interpret parentheses and quotes in dnl comments, and
dnl  having a whole file as a macro argument would overflow the string space
dnl  on BSD m4.

AC_DEFUN([GMP_FINISH],
[AC_REQUIRE([GMP_INIT])
echo "creating $gmp_configm4"
echo ["d""nl $gmp_configm4.  Generated automatically by configure."] > $gmp_configm4
if test -f $gmp_tmpconfigm4; then
  echo ["changequote(<,>)"] >> $gmp_configm4
  echo ["ifdef(<__CONFIG_M4_INCLUDED__>,,<"] >> $gmp_configm4
  cat $gmp_tmpconfigm4 >> $gmp_configm4
  echo [">)"] >> $gmp_configm4
  echo ["changequote(\`,')"] >> $gmp_configm4
  rm $gmp_tmpconfigm4
fi
echo ["ifdef(\`__CONFIG_M4_INCLUDED__',,\`"] >> $gmp_configm4
if test -f $gmp_tmpconfigm4i; then
  cat $gmp_tmpconfigm4i >> $gmp_configm4
  rm $gmp_tmpconfigm4i
fi
if test -f $gmp_tmpconfigm4p; then
  cat $gmp_tmpconfigm4p >> $gmp_configm4
  rm $gmp_tmpconfigm4p
fi
echo ["')"] >> $gmp_configm4
echo ["define(\`__CONFIG_M4_INCLUDED__')"] >> $gmp_configm4
])


dnl  GMP_INCLUDE_MPN(FILE)
dnl  ---------------------
dnl  Add an include_mpn(`FILE') to config.m4.  FILE should be a path
dnl  relative to the mpn source directory, for example
dnl
dnl      GMP_INCLUDE_MPN(`x86/x86-defs.m4')
dnl

AC_DEFUN([GMP_INCLUDE_MPN],
[AC_REQUIRE([GMP_INIT])
echo ["include][(CONFIG_TOP_SRCDIR\`/$1')"] >>$gmp_tmpconfigm4i
])


dnl  GMP_PROG_M4
dnl  -----------
dnl  Find a working m4, either in $PATH or likely locations, and setup $M4
dnl  and an AC_SUBST accordingly.  If $M4 is already set then it's a user
dnl  choice and is accepted with no checks.  GMP_PROG_M4 is like
dnl  AC_PATH_PROG or AC_CHECK_PROG, but tests each m4 found to see if it's
dnl  good enough.
dnl
dnl  See mpn/asm-defs.m4 for details on the known bad m4s.

AC_DEFUN([GMP_PROG_M4],
[AC_ARG_VAR(M4,[m4 macro processor])
AC_CACHE_CHECK([for suitable m4],
                gmp_cv_prog_m4,
[if test -n "$M4"; then
  gmp_cv_prog_m4="$M4"
else
  cat >conftest.m4 <<\EOF
dnl  Must protect this against being expanded during autoconf m4!
dnl  Dont put "dnl"s in this as autoconf will flag an error for unexpanded
dnl  macros.
[define(dollarhash,``$][#'')ifelse(dollarhash(x),1,`define(t1,Y)',
``bad: $][# not supported (SunOS /usr/bin/m4)
'')ifelse(eval(89),89,`define(t2,Y)',
`bad: eval() doesnt support 8 or 9 in a constant (OpenBSD 2.6 m4)
')ifelse(eval(9,9),10,`define(t3,Y)',
`bad: eval() doesnt support radix in eval (FreeBSD 8.x,9.0,9.1,9.2 m4)
')ifelse(t1`'t2`'t3,YYY,`good
')]
EOF
dnl ' <- balance the quotes for emacs sh-mode
  echo "trying m4" >&AS_MESSAGE_LOG_FD
  gmp_tmp_val=`(m4 conftest.m4) 2>&AS_MESSAGE_LOG_FD`
  echo "$gmp_tmp_val" >&AS_MESSAGE_LOG_FD
  if test "$gmp_tmp_val" = good; then
    gmp_cv_prog_m4="m4"
  else
    IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS=":"
dnl $ac_dummy forces splitting on constant user-supplied paths.
dnl POSIX.2 word splitting is done only on the output of word expansions,
dnl not every word.  This closes a longstanding sh security hole.
    ac_dummy="$PATH:/usr/5bin"
    for ac_dir in $ac_dummy; do
      test -z "$ac_dir" && ac_dir=.
      echo "trying $ac_dir/m4" >&AS_MESSAGE_LOG_FD
      gmp_tmp_val=`($ac_dir/m4 conftest.m4) 2>&AS_MESSAGE_LOG_FD`
      echo "$gmp_tmp_val" >&AS_MESSAGE_LOG_FD
      if test "$gmp_tmp_val" = good; then
        gmp_cv_prog_m4="$ac_dir/m4"
        break
      fi
    done
    IFS="$ac_save_ifs"
    if test -z "$gmp_cv_prog_m4"; then
      AC_MSG_ERROR([No usable m4 in \$PATH or /usr/5bin (see config.log for reasons).])
    fi
  fi
  rm -f conftest.m4
fi])
M4="$gmp_cv_prog_m4"
AC_SUBST(M4)
])

dnl  GMP_PATH_NM
dnl  -----------
dnl  GMP additions to libtool LT_PATH_NM.
dnl
dnl  Note that if LT_PATH_NM can't find a working nm it still leaves
dnl  $NM set to "nm", so $NM can't be assumed to actually work.
dnl
dnl  A user-selected $NM is always left unchanged.  LT_PATH_NM is still run
dnl  to get the "checking" message printed though.
dnl
dnl  Perhaps it'd be worthwhile checking that nm works, by running it on an
dnl  actual object file.  For instance on sparcv9 solaris old versions of
dnl  GNU nm don't recognise 64-bit objects.  Checking would give a better
dnl  error message than just a failure in later tests like GMP_ASM_W32 etc.
dnl
dnl  On the other hand it's not really normal autoconf practice to take too
dnl  much trouble over detecting a broken set of tools.  And libtool doesn't
dnl  do anything at all for say ranlib or strip.  So for now we're inclined
dnl  to just demand that the user provides a coherent environment.

AC_DEFUN([GMP_PATH_NM],
[
gmp_user_NM=$NM
dnl Is already called in LT_INIT
dnl LT_PATH_NM

# FIXME: When cross compiling (ie. $ac_tool_prefix not empty), libtool
# defaults to plain "nm" if a "${ac_tool_prefix}nm" is not found.  In this
# case run it again to try the native "nm", firstly so that likely locations
# are searched, secondly so that -B or -p are added if necessary for BSD
# format.  This is necessary for instance on OSF with "./configure
# --build=alphaev5-dec-osf --host=alphaev6-dec-osf".
#
if test -z "$gmp_user_NM" && test -n "$ac_tool_prefix" && test "$NM" = nm; then
  $as_unset lt_cv_path_NM
  gmp_save_ac_tool_prefix=$ac_tool_prefix
  ac_tool_prefix=
  NM=
  LT_PATH_NM
  ac_tool_prefix=$gmp_save_ac_tool_prefix
fi

if test -z "$gmp_user_NM"; then
                        eval nmflags=\"\$nm${abi1}_flags\"
  test -n "$nmflags" || eval nmflags=\"\$nm${abi2}_flags\"
  if test -n "$nmflags"; then
    AC_MSG_CHECKING([for extra nm flags])
    NM="$NM $nmflags"
    AC_MSG_RESULT([$nmflags])
  fi
fi
])


dnl  GMP_DEFINE(MACRO, DEFINITION [, LOCATION])
dnl  ------------------------------------------
dnl  Define M4 macro MACRO as DEFINITION in temporary file.
dnl
dnl  If LOCATION is `POST', the definition will appear after any include()
dnl  directives inserted by GMP_INCLUDE.  Mind the quoting!  No shell
dnl  variables will get expanded.  Don't forget to invoke GMP_FINISH to
dnl  create file config.m4.  config.m4 uses `<' and '>' as quote characters
dnl  for all defines.

AC_DEFUN([GMP_DEFINE],
[AC_REQUIRE([GMP_INIT])
echo ['define(<$1>, <$2>)'] >>ifelse([$3], [POST],
                              $gmp_tmpconfigm4p, $gmp_tmpconfigm4)
])


dnl  GMP_DEFINE_RAW(STRING [, LOCATION])
dnl  ------------------------------------
dnl  Put STRING into config.m4 file.
dnl
dnl  If LOCATION is `POST', the definition will appear after any include()
dnl  directives inserted by GMP_INCLUDE.  Don't forget to invoke GMP_FINISH
dnl  to create file config.m4.

AC_DEFUN([GMP_DEFINE_RAW],
[AC_REQUIRE([GMP_INIT])
echo [$1] >> ifelse([$2], [POST], $gmp_tmpconfigm4p, $gmp_tmpconfigm4)
])


dnl  GMP_TRY_ASSEMBLE(asm-code,[action-success][,action-fail])
dnl  ----------------------------------------------------------
dnl  Attempt to assemble the given code.
dnl  Do "action-success" if this succeeds, "action-fail" if not.
dnl
dnl  conftest.o and conftest.out are available for inspection in
dnl  "action-success".  If either action does a "break" out of a loop then
dnl  an explicit "rm -f conftest*" will be necessary.
dnl
dnl  This is not unlike AC_TRY_COMPILE, but there's no default includes or
dnl  anything in "asm-code", everything wanted must be given explicitly.

AC_DEFUN([GMP_TRY_ASSEMBLE],
[cat >conftest.s <<EOF
[$1]
EOF
gmp_assemble="$CC -c $CFLAGS $CPPFLAGS conftest.s >conftest.out 2>&1"
if AC_TRY_EVAL(gmp_assemble); then
  cat conftest.out >&AS_MESSAGE_LOG_FD
  ifelse([$2],,:,[$2])
else
  cat conftest.out >&AS_MESSAGE_LOG_FD
  echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
  cat conftest.s >&AS_MESSAGE_LOG_FD
  ifelse([$3],,:,[$3])
fi
rm -f conftest*
])


dnl  CL_ASM_NOEXECSTACK
dnl  -------------------
dnl
dnl  Checks whether the stack can be marked nonexecutable by passing an option
dnl  to the C-compiler when acting on .s files. Appends that option to ASMFLAGS.
dnl  This macro is adapted from one found in GLIBC-2.3.5.
dnl
dnl  FIXME: This test looks broken. It tests that a file with
dnl  .note.GNU-stack... can be compiled/assembled with -Wa,--noexecstack.  It
dnl  does not determine if that command-line option has any effect on general
dnl  asm code.
AC_DEFUN([CL_ASM_NOEXECSTACK],
[AC_REQUIRE([AC_PROG_CC])
AC_CACHE_CHECK([whether assembler supports --noexecstack option],
                cl_cv_asm_noexecstack,
[dnl
  cat > conftest.c <<EOF
void foo() {}
EOF
  if AC_TRY_COMMAND([${CC} $CFLAGS $CPPFLAGS
                     -S -o conftest.s conftest.c >/dev/null]) \
     && grep .note.GNU-stack conftest.s >/dev/null \
     && AC_TRY_COMMAND([${CC} $CFLAGS $CPPFLAGS -Wa,--noexecstack
                       -c -o conftest.o conftest.s >/dev/null])
  then
    cl_cv_asm_noexecstack=yes
  else
    cl_cv_asm_noexecstack=no
  fi
  rm -f conftest*])
  if test "$cl_cv_asm_noexecstack" = yes; then
    ASMFLAGS="$ASMFLAGS -Wa,--noexecstack"
  fi
  AC_SUBST(ASMFLAGS)
])


dnl  GMP_ASM_LABEL_SUFFIX
dnl  --------------------
dnl  : - is usual.
dnl  empty - hppa on HP-UX doesn't use a :, just the label name
dnl
dnl  Note that it's necessary to test the empty case first, since HP "as"
dnl  will accept "somelabel:", and take it to mean a label with a name that
dnl  happens to end in a colon.

AC_DEFUN([GMP_ASM_LABEL_SUFFIX],
[AC_REQUIRE([GMP_ASM_TEXT])
AC_CACHE_CHECK([for assembler label suffix],
                gmp_cv_asm_label_suffix,
[gmp_cv_asm_label_suffix=unknown
for i in "" ":"; do
  echo "trying $i" >&AS_MESSAGE_LOG_FD
  GMP_TRY_ASSEMBLE(
[	$gmp_cv_asm_text
somelabel$i],
    [gmp_cv_asm_label_suffix=$i
     rm -f conftest*
     break],
    [cat conftest.out >&AS_MESSAGE_LOG_FD])
done
if test "$gmp_cv_asm_label_suffix" = "unknown"; then
  AC_MSG_ERROR([Cannot determine label suffix])
fi
])
echo ["define(<LABEL_SUFFIX>, <$gmp_cv_asm_label_suffix>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_UNDERSCORE
dnl  ------------------
dnl  Determine whether global symbols need to be prefixed with an underscore.
dnl  The output from "nm" is grepped to see what a typical symbol looks like.
dnl
dnl  This test used to grep the .o file directly, but that failed with greps
dnl  that don't like binary files (eg. SunOS 4).
dnl
dnl  This test also used to construct an assembler file with and without an
dnl  underscore and try to link that to a C file, to see which worked.
dnl  Although that's what will happen in the real build we don't really want
dnl  to depend on creating asm files within configure for every possible CPU
dnl  (or at least we don't want to do that more than we have to).
dnl
dnl  The fallback on no underscore is based on the assumption that the world
dnl  is moving towards non-underscore systems.  There should actually be no
dnl  good reason for nm to fail though.

AC_DEFUN([GMP_ASM_UNDERSCORE],
[AC_REQUIRE([GMP_PATH_NM])
AC_CACHE_CHECK([if globals are prefixed by underscore],
               gmp_cv_asm_underscore,
[gmp_cv_asm_underscore="unknown"
cat >conftest.c <<EOF
int gurkmacka;
EOF
gmp_compile="$CC $CFLAGS $CPPFLAGS -c conftest.c >&AS_MESSAGE_LOG_FD"
if AC_TRY_EVAL(gmp_compile); then
  $NM conftest.$OBJEXT >conftest.out
  if grep "[[ 	]]_gurkmacka" conftest.out >/dev/null; then
    gmp_cv_asm_underscore=yes
  elif grep "[[ 	]]gurkmacka" conftest.out >/dev/null; then
    gmp_cv_asm_underscore=no
  else
    echo "configure: $NM doesn't have gurkmacka:" >&AS_MESSAGE_LOG_FD
    cat conftest.out >&AS_MESSAGE_LOG_FD
  fi
else
  echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
  cat conftest.c >&AS_MESSAGE_LOG_FD
fi
rm -f conftest*
])
case $gmp_cv_asm_underscore in
  yes)
    GMP_DEFINE(GSYM_PREFIX, [_]) ;;
  no)
    GMP_DEFINE(GSYM_PREFIX, []) ;;
  *)
    AC_MSG_WARN([+----------------------------------------------------------])
    AC_MSG_WARN([| Cannot determine global symbol prefix.])
    AC_MSG_WARN([| $NM output doesn't contain a global data symbol.])
    AC_MSG_WARN([| Will proceed with no underscore.])
    AC_MSG_WARN([| If this is wrong then you'll get link errors referring])
    AC_MSG_WARN([| to ___gmpn_add_n (note three underscores).])
    AC_MSG_WARN([| In this case do a fresh build with an override,])
    AC_MSG_WARN([|     ./configure gmp_cv_asm_underscore=yes])
    AC_MSG_WARN([+----------------------------------------------------------])
    GMP_DEFINE(GSYM_PREFIX, [])
    ;;
esac
])


dnl  GMP_ASM_ALIGN_LOG
dnl  -----------------
dnl  Is parameter to `.align' logarithmic?

AC_DEFUN([GMP_ASM_ALIGN_LOG],
[AC_REQUIRE([GMP_ASM_GLOBL])
AC_REQUIRE([GMP_ASM_BYTE])
AC_REQUIRE([GMP_ASM_DATA])
AC_REQUIRE([GMP_ASM_LABEL_SUFFIX])
AC_REQUIRE([GMP_PATH_NM])
AC_CACHE_CHECK([if .align assembly directive is logarithmic],
               gmp_cv_asm_align_log,
[GMP_TRY_ASSEMBLE(
[	$gmp_cv_asm_data
	.align  4
	$gmp_cv_asm_globl	foo
	$gmp_cv_asm_byte	1
	.align	4
foo$gmp_cv_asm_label_suffix
	$gmp_cv_asm_byte	2],
  [gmp_tmp_val=[`$NM conftest.$OBJEXT | grep foo | \
     sed -e 's;[[][0-9][]]\(.*\);\1;' -e 's;[^1-9]*\([0-9]*\).*;\1;'`]
  if test "$gmp_tmp_val" = "10" || test "$gmp_tmp_val" = "16"; then
    gmp_cv_asm_align_log=yes
  else
    gmp_cv_asm_align_log=no
  fi],
  [AC_MSG_ERROR([cannot assemble alignment test])])])

GMP_DEFINE_RAW(["define(<ALIGN_LOGARITHMIC>,<$gmp_cv_asm_align_log>)"])
])


dnl  GMP_ASM_ALIGN_FILL_0x90
dnl  -----------------------
dnl  Determine whether a ",0x90" suffix works on a .align directive.
dnl  This is only meant for use on x86, 0x90 being a "nop".
dnl
dnl  Old gas, eg. 1.92.3
dnl       Needs ",0x90" or else the fill is 0x00, which can't be executed
dnl       across.
dnl
dnl  New gas, eg. 2.91
dnl       Generates multi-byte nop fills even when ",0x90" is given.
dnl
dnl  Solaris 2.6 as
dnl       ",0x90" is not allowed, causes a fatal error.
dnl
dnl  Solaris 2.8 as
dnl       ",0x90" does nothing, generates a warning that it's being ignored.
dnl
dnl  SCO OpenServer 5 as
dnl       Second parameter is max bytes to fill, not a fill pattern.
dnl       ",0x90" is an error due to being bigger than the first parameter.
dnl       Multi-byte nop fills are generated in text segments.
dnl
dnl  Note that both solaris "as"s only care about ",0x90" if they actually
dnl  have to use it to fill something, hence the .byte in the test.  It's
dnl  the second .align which provokes the error or warning.
dnl
dnl  The warning from solaris 2.8 is suppressed to stop anyone worrying that
dnl  something might be wrong.

AC_DEFUN([GMP_ASM_ALIGN_FILL_0x90],
[AC_REQUIRE([GMP_ASM_TEXT])
AC_CACHE_CHECK([if the .align directive accepts an 0x90 fill in .text],
               gmp_cv_asm_align_fill_0x90,
[GMP_TRY_ASSEMBLE(
[	$gmp_cv_asm_text
	.align  4, 0x90
	.byte   0
	.align  4, 0x90],
[if grep "Warning: Fill parameter ignored for executable section" conftest.out >/dev/null; then
  echo "Suppressing this warning by omitting 0x90" 1>&AS_MESSAGE_LOG_FD
  gmp_cv_asm_align_fill_0x90=no
else
  gmp_cv_asm_align_fill_0x90=yes
fi],
[gmp_cv_asm_align_fill_0x90=no])])

GMP_DEFINE_RAW(["define(<ALIGN_FILL_0x90>,<$gmp_cv_asm_align_fill_0x90>)"])
])


dnl  GMP_ASM_BYTE
dnl  ------------
dnl  .byte - is usual.
dnl  data1 - required by ia64 (on hpux at least).
dnl
dnl  This macro is just to support other configure tests, not any actual asm
dnl  code.

AC_DEFUN([GMP_ASM_BYTE],
[AC_REQUIRE([GMP_ASM_TEXT])
AC_REQUIRE([GMP_ASM_LABEL_SUFFIX])
AC_CACHE_CHECK([for assembler byte directive],
                gmp_cv_asm_byte,
[for i in .byte data1; do
  echo "trying $i" >&AS_MESSAGE_LOG_FD
  GMP_TRY_ASSEMBLE(
[	$gmp_cv_asm_data
	$i	0
],
    [gmp_cv_asm_byte=$i
     rm -f conftest*
     break],
    [cat conftest.out >&AS_MESSAGE_LOG_FD])
done
if test -z "$gmp_cv_asm_byte"; then
  AC_MSG_ERROR([Cannot determine how to emit a data byte])
fi
])
])


dnl  GMP_ASM_TEXT
dnl  ------------
dnl  .text - is usual.
dnl  .code - is needed by the hppa on HP-UX (but ia64 HP-UX uses .text)
dnl  .csect .text[PR] - is for AIX.

AC_DEFUN([GMP_ASM_TEXT],
[AC_CACHE_CHECK([how to switch to text section],
                gmp_cv_asm_text,
[for i in ".text" ".code" [".csect .text[PR]"]; do
  echo "trying $i" >&AS_MESSAGE_LOG_FD
  GMP_TRY_ASSEMBLE([	$i],
    [gmp_cv_asm_text=$i
     rm -f conftest*
     break])
done
if test -z "$gmp_cv_asm_text"; then
  AC_MSG_ERROR([Cannot determine text section directive])
fi
])
echo ["define(<TEXT>, <$gmp_cv_asm_text>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_DATA
dnl  ------------
dnl  Can we say `.data'?

AC_DEFUN([GMP_ASM_DATA],
[AC_CACHE_CHECK([how to switch to data section],
                gmp_cv_asm_data,
[case $host in
  *-*-aix*) gmp_cv_asm_data=[".csect .data[RW]"] ;;
  *)        gmp_cv_asm_data=".data" ;;
esac
])
echo ["define(<DATA>, <$gmp_cv_asm_data>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_RODATA
dnl  --------------
dnl  Find out how to switch to the read-only data section.
dnl
dnl  The compiler output is grepped for the right directive.  It's not
dnl  considered wise to just probe for ".section .rodata" or whatever works,
dnl  since arbitrary section names might be accepted, but not necessarily do
dnl  the right thing when they get to the linker.
dnl
dnl  Only a few asm files use RODATA, so this code is perhaps a bit
dnl  excessive right now, but should find more uses in the future.
dnl
dnl  FIXME: gcc on aix generates something like ".csect _foo.ro_c[RO],3"
dnl  where foo is the object file.  Might need to check for that if we use
dnl  RODATA there.

AC_DEFUN([GMP_ASM_RODATA],
[AC_REQUIRE([GMP_ASM_TEXT])
AC_REQUIRE([GMP_ASM_DATA])
AC_REQUIRE([GMP_ASM_LABEL_SUFFIX])
AC_REQUIRE([GMP_ASM_UNDERSCORE])
AC_CACHE_CHECK([how to switch to read-only data section],
               gmp_cv_asm_rodata,
[
dnl Default to DATA on CPUs with split code/data caching, and TEXT
dnl elsewhere.  i386 means generic x86, so use DATA on it.
case $host in
X86_64_PATTERN)
  gmp_cv_asm_rodata="$gmp_cv_asm_data" ;;
*)
  gmp_cv_asm_rodata="$gmp_cv_asm_text" ;;
esac

cat >conftest.c <<EOF
extern const int foo[[]];		/* Suppresses C++'s suppression of foo */
const int foo[[]] = {1,2,3};
EOF
echo "Test program:" >&AS_MESSAGE_LOG_FD
cat conftest.c >&AS_MESSAGE_LOG_FD
gmp_compile="$CC $CFLAGS $CPPFLAGS -S conftest.c >&AS_MESSAGE_LOG_FD"
if AC_TRY_EVAL(gmp_compile); then
  echo "Compiler output:" >&AS_MESSAGE_LOG_FD
  cat conftest.s >&AS_MESSAGE_LOG_FD
  if test $gmp_cv_asm_underscore = yes; then
    tmp_gsym_prefix=_
  else
    tmp_gsym_prefix=
  fi
  # must see our label
  if grep "^${tmp_gsym_prefix}foo$gmp_cv_asm_label_suffix" conftest.s >/dev/null 2>&AS_MESSAGE_LOG_FD; then
    # take the last directive before our label (hence skipping segments
    # getting debugging info etc)
    tmp_match=`sed -n ["/^${tmp_gsym_prefix}foo$gmp_cv_asm_label_suffix/q
                        /^[. 	]*data/p
                        /^[. 	]*rdata/p
                        /^[. 	]*text/p
                        /^[. 	]*section/p
                        /^[. 	]*csect/p
                        /^[. 	]*CSECT/p"] conftest.s | sed -n '$p'`
    echo "Match: $tmp_match" >&AS_MESSAGE_LOG_FD
    if test -n "$tmp_match"; then
      gmp_cv_asm_rodata=$tmp_match
    fi
  else
    echo "Couldn't find label: ^${tmp_gsym_prefix}foo$gmp_cv_asm_label_suffix" >&AS_MESSAGE_LOG_FD
  fi
fi
rm -f conftest*
])
echo ["define(<RODATA>, <$gmp_cv_asm_rodata>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_GLOBL
dnl  -------------
dnl  The assembler directive to mark a label as a global symbol.
dnl
dnl  ia64 - .global is standard, according to the Intel documentation.
dnl
dnl  hppa - ".export foo,entry" is demanded by HP hppa "as".  ".global" is a
dnl      kind of import.
dnl
dnl  other - .globl is usual.
dnl
dnl  "gas" tends to accept .globl everywhere, in addition to .export or
dnl  .global or whatever the system assembler demands.

AC_DEFUN([GMP_ASM_GLOBL],
[AC_REQUIRE([GMP_ASM_TEXT])
AC_CACHE_CHECK([for assembler global directive],
                gmp_cv_asm_globl,
[case $host in
  dnl We do not support HPPA or IA64 assembly
  dnl hppa*-*-*)     gmp_cv_asm_globl=.export ;;
  dnl IA64_PATTERN)  gmp_cv_asm_globl=.global ;;
  *)             gmp_cv_asm_globl=.globl  ;;
esac
])
echo ["define(<GLOBL>, <$gmp_cv_asm_globl>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_GLOBL_ATTR
dnl  ------------------
dnl  Do we need something after `GLOBL symbol'?

AC_DEFUN([GMP_ASM_GLOBL_ATTR],
[AC_REQUIRE([GMP_ASM_GLOBL])
AC_CACHE_CHECK([for assembler global directive attribute],
                gmp_cv_asm_globl_attr,
[case $gmp_cv_asm_globl in
  .export) gmp_cv_asm_globl_attr=",entry" ;;
  *)       gmp_cv_asm_globl_attr="" ;;
esac
])
echo ["define(<GLOBL_ATTR>, <$gmp_cv_asm_globl_attr>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_TYPE
dnl  ------------
dnl  Can we say ".type", and how?
dnl
dnl  For i386 GNU/Linux ELF systems, and very likely other ELF systems,
dnl  .type and .size are important on functions in shared libraries.  If
dnl  .type is omitted and the mainline program references that function then
dnl  the code will be copied down to the mainline at load time like a piece
dnl  of data.  If .size is wrong or missing (it defaults to 4 bytes or some
dnl  such) then incorrect bytes will be copied and a segv is the most likely
dnl  result.  In any case such copying is not what's wanted, a .type
dnl  directive will ensure a PLT entry is used.
dnl
dnl  In GMP the assembler functions are normally only used from within the
dnl  library (since most programs are not interested in the low level
dnl  routines), and in those circumstances a missing .type isn't fatal,
dnl  letting the problem go unnoticed.  tests/mpn/t-asmtype.c aims to check
dnl  for it.

AC_DEFUN([GMP_ASM_TYPE],
[AC_CACHE_CHECK([for assembler .type directive],
                gmp_cv_asm_type,
[gmp_cv_asm_type=
for gmp_tmp_prefix in @ \# %; do
  GMP_TRY_ASSEMBLE([	.type	sym,${gmp_tmp_prefix}function],
    [if grep "\.type pseudo-op used outside of \.def/\.endef ignored" conftest.out >/dev/null; then : ;
    else
      gmp_cv_asm_type=".type	\$][1,${gmp_tmp_prefix}\$][2"
      break
    fi])
done
rm -f conftest*
])
echo ["define(<TYPE>, <$gmp_cv_asm_type>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_SIZE
dnl  ------------
dnl  Can we say `.size'?

AC_DEFUN([GMP_ASM_SIZE],
[AC_CACHE_CHECK([for assembler .size directive],
                gmp_cv_asm_size,
[gmp_cv_asm_size=
GMP_TRY_ASSEMBLE([	.size	sym,1],
  [if grep "\.size pseudo-op used outside of \.def/\.endef ignored" conftest.out >/dev/null; then : ;
  else
    gmp_cv_asm_size=".size	\$][1,\$][2"
  fi])
])
echo ["define(<SIZE>, <$gmp_cv_asm_size>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_COFF_TYPE
dnl  -----------------
dnl  Determine whether the assembler supports COFF type information.
dnl
dnl  Currently this is only needed for mingw (and cygwin perhaps) and so is
dnl  run only on the x86s, but it ought to work anywhere.
dnl
dnl  On MINGW, recent versions of the linker have an automatic import scheme
dnl  for data in a DLL which is referenced by a mainline but without
dnl  __declspec (__dllimport__) on the prototype.  It seems functions
dnl  without type information are treated as data, or something, and calls
dnl  to them from the mainline will crash.  gcc puts type information on the
dnl  C functions it generates, we need to do the same for assembler
dnl  functions.
dnl
dnl  This applies only to functions without __declspec(__dllimport__),
dnl  ie. without __GMP_DECLSPEC in the case of libgmp, so it also works just
dnl  to ensure all assembler functions used from outside libgmp have
dnl  __GMP_DECLSPEC on their prototypes.  But this isn't an ideal situation,
dnl  since we don't want perfectly valid calls going wrong just because
dnl  there wasn't a prototype in scope.
dnl
dnl  When an auto-import takes place, the following warning is given by the
dnl  linker.  This shouldn't be seen for any functions.
dnl
dnl      Info: resolving _foo by linking to __imp__foo (auto-import)
dnl
dnl
dnl  COFF type directives look like the following
dnl
dnl      .def    _foo
dnl      .scl    2
dnl      .type   32
dnl      .endef
dnl
dnl  _foo is the symbol with GSYM_PREFIX (_).  .scl is the storage class, 2
dnl  for external, 3 for static.  .type is the object type, 32 for a
dnl  function.
dnl
dnl  On an ELF system, this is (correctly) rejected due to .def, .endef and
dnl  .scl being invalid, and .type not having enough arguments.

AC_DEFUN([GMP_ASM_COFF_TYPE],
[AC_REQUIRE([GMP_ASM_TEXT])
AC_REQUIRE([GMP_ASM_GLOBL])
AC_REQUIRE([GMP_ASM_GLOBL_ATTR])
AC_REQUIRE([GMP_ASM_LABEL_SUFFIX])
AC_REQUIRE([GMP_ASM_UNDERSCORE])
AC_CACHE_CHECK([for assembler COFF type directives],
		gmp_cv_asm_x86_coff_type,
[GMP_TRY_ASSEMBLE(
[	$gmp_cv_asm_text
	$gmp_cv_asm_globl ${tmp_gsym_prefix}foo$gmp_cv_asm_globl_attr
	.def	${tmp_gsym_prefix}foo
	.scl	2
	.type	32
	.endef
${tmp_gsym_prefix}foo$gmp_cv_asm_label_suffix
],
  [gmp_cv_asm_x86_coff_type=yes],
  [gmp_cv_asm_x86_coff_type=no])
])
echo ["define(<HAVE_COFF_TYPE>, <$gmp_cv_asm_x86_coff_type>)"] >> $gmp_tmpconfigm4
])


dnl  GMP_ASM_LSYM_PREFIX
dnl  -------------------
dnl  What is the prefix for a local label?
dnl
dnl  The prefixes tested are,
dnl
dnl      L  - usual for underscore systems
dnl      .L - usual for non-underscore systems
dnl      $  - alpha (gas and OSF system assembler)
dnl      L$ - hppa (gas and HP-UX system assembler)
dnl
dnl  The default is "L" if the tests fail for any reason.  There's a good
dnl  chance this will be adequate, since on most systems labels are local
dnl  anyway unless given a ".globl", and an "L" will avoid clashes with
dnl  other identifiers.
dnl
dnl  For gas, ".L" is normally purely local to the assembler, it doesn't get
dnl  put into the object file at all.  This style is preferred, to keep the
dnl  object files nice and clean.
dnl
dnl  BSD format nm produces a line like
dnl
dnl      00000000 t Lgurkmacka
dnl
dnl  The symbol code is normally "t" for text, but any lower case letter
dnl  indicates a local definition.
dnl
dnl  Code "n" is for a debugging symbol, OSF "nm -B" gives that as an upper
dnl  case "N" for a local.
dnl
dnl  HP-UX nm prints an error message (though seems to give a 0 exit) if
dnl  there's no symbols at all in an object file, hence the use of "dummy".

AC_DEFUN([GMP_ASM_LSYM_PREFIX],
[AC_REQUIRE([GMP_ASM_LABEL_SUFFIX])
AC_REQUIRE([GMP_ASM_TEXT])
AC_REQUIRE([GMP_PATH_NM])
AC_CACHE_CHECK([for assembler local label prefix],
               gmp_cv_asm_lsym_prefix,
[gmp_tmp_pre_appears=yes
for gmp_tmp_pre in L .L $L $ L$; do
  echo "Trying $gmp_tmp_pre" >&AS_MESSAGE_LOG_FD
  GMP_TRY_ASSEMBLE(
[	$gmp_cv_asm_text
dummy${gmp_cv_asm_label_suffix}
${gmp_tmp_pre}gurkmacka${gmp_cv_asm_label_suffix}],
  [if $NM conftest.$OBJEXT >conftest.nm 2>&AS_MESSAGE_LOG_FD; then : ; else
    cat conftest.nm >&AS_MESSAGE_LOG_FD
    AC_MSG_WARN(["$NM" failure])
    break
  fi
  cat conftest.nm >&AS_MESSAGE_LOG_FD
  if grep gurkmacka conftest.nm >/dev/null; then : ; else
    # no mention of the symbol, this is good
    echo "$gmp_tmp_pre label doesn't appear in object file at all (good)" >&AS_MESSAGE_LOG_FD
    gmp_cv_asm_lsym_prefix="$gmp_tmp_pre"
    gmp_tmp_pre_appears=no
    break
  fi
  if grep [' [a-zN] .*gurkmacka'] conftest.nm >/dev/null; then
    # symbol mentioned as a local, use this if nothing better
    echo "$gmp_tmp_pre label is local but still in object file" >&AS_MESSAGE_LOG_FD
    if test -z "$gmp_cv_asm_lsym_prefix"; then
      gmp_cv_asm_lsym_prefix="$gmp_tmp_pre"
    fi
  else
    echo "$gmp_tmp_pre label is something unknown" >&AS_MESSAGE_LOG_FD
  fi
  ])
done
rm -f conftest*
if test -z "$gmp_cv_asm_lsym_prefix"; then
  gmp_cv_asm_lsym_prefix=L
  AC_MSG_WARN([cannot determine local label, using default $gmp_cv_asm_lsym_prefix])
fi
# for development purposes, note whether we got a purely temporary local label
echo "Local label appears in object files: $gmp_tmp_pre_appears" >&AS_MESSAGE_LOG_FD
])
echo ["define(<LSYM_PREFIX>, <${gmp_cv_asm_lsym_prefix}>)"] >> $gmp_tmpconfigm4
AC_DEFINE_UNQUOTED(LSYM_PREFIX, "$gmp_cv_asm_lsym_prefix",
                   [Assembler local label prefix])
])
