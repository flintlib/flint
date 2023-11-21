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

# AX_CHECK_LONGLONG_HEADER(FLAG, [ACTION-SUCCESS], [ACTION-FAILURE], [EXTRA-PROLOGUE])
AC_DEFUN([AX_CHECK_LONGLONG_HEADER],
[AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_$1])dnl
AC_CACHE_CHECK([whether system can compile with $1], CACHEVAR, [
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
    [[#include <gmp.h>
      typedef mp_limb_t ulong;
      #define FLINT_ASSERT(x)
      #define FLINT_DLL
      #define FLINT_BITS GMP_LIMB_BITS
      #include "$1"
      $4]],
    [[ulong s3, s2, s1, s0;
      ulong a3, a2, a1, a0;
      ulong b3, b2, b1, b0;

      a0 = 19827;
      a1 = 1872;
      a2 = 338237;
      a3 = 98080;
      b0 = 1798291;
      b1 = 719271;
      b2 = 891;
      b3 = 9112;

      s0 = flint_clz(a0);
      s1 = flint_ctz(a0);
      add_ssaaaa(s1, s0, a1, a0, b1, b0);
      add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0);
      add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0);
      sub_ddmmss(s1, s0, a1, a0, b1, b0);
      sub_dddmmmsss(s2, s1, s0, a2, a1, a0, b2, b1, b0);
      umul_ppmm(s1, s0, a0, b0);
      smul_ppmm(s1, s0, a0, b0);
      udiv_qrnnd(s1, s0, a1, a0, b0);
      sdiv_qrnnd(s1, s0, a1, a0, b0);
      byte_swap(s2);
      /* udiv_qrnnd_preinv is the same for every header */]]
    )],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])])
AS_VAR_IF(CACHEVAR,yes,
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_LONGLONG_HEADER
