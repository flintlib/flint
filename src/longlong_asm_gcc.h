/*
    Copyright 1991-1994, 1996, 1997, 1999-2005, 2007-2009, 2011-2020 Free
    Software Foundation, Inc.

    Copyright (C) 2023, 2024 Albin Ahlbäck

    This file is part of FLINT.

    Contains code from GNU MP Library.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LONGLONG_ASM_H
#define LONGLONG_ASM_H

/* Machine specific operations */
#if defined (__amd64__) || (FLINT_BITS == 32 && (defined (__i386__) || defined (__i486__)))

# if FLINT_BITS == 64 && defined (__amd64__)
#  define _ASM_ADD "addq"
#  define _ASM_ADC "adcq"
#  define _ASM_SUB "subq"
#  define _ASM_SBB "sbbq"
#  define _ASM_MUL "mulq"
#  define _ASM_IMUL "imulq"
#  define _ASM_PRE "q"
#  define _ASM_RME "rme"
# else
#  define _ASM_ADD "addl"
#  define _ASM_ADC "adcl"
#  define _ASM_SUB "subl"
#  define _ASM_SBB "sbbl"
#  define _ASM_MUL "mull"
#  define _ASM_IMUL "imull"
#  define _ASM_PRE "k"
#  define _ASM_RME "g"
# endif

# define add_ssaaaa(s1, s0, a1, a0, b1, b0) \
  __asm__(_ASM_ADD " %5,%" _ASM_PRE "1\n" \
     "\t" _ASM_ADC " %3,%" _ASM_PRE "0" \
    : "=r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "%1" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  __asm__(_ASM_ADD " %8,%" _ASM_PRE "2\n" \
     "\t" _ASM_ADC " %6,%" _ASM_PRE "1\n" \
     "\t" _ASM_ADC " %4,%" _ASM_PRE "0" \
    : "=r" (s2), "=&r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a2)), _ASM_RME ((ulong)(b2)), \
      "1" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "2" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  __asm__(_ASM_ADD " %11,%" _ASM_PRE "3\n" \
     "\t" _ASM_ADC " %9,%" _ASM_PRE "2\n" \
     "\t" _ASM_ADC " %7,%" _ASM_PRE "1\n" \
     "\t" _ASM_ADC " %5,%" _ASM_PRE "0" \
    : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a3)), _ASM_RME ((ulong)(b3)), \
      "1" ((ulong)(a2)), _ASM_RME ((ulong)(b2)), \
      "2" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "3" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0) \
  __asm__(_ASM_ADD " %14,%" _ASM_PRE "4\n" \
     "\t" _ASM_ADC " %12,%" _ASM_PRE "3\n" \
     "\t" _ASM_ADC " %10,%" _ASM_PRE "2\n" \
     "\t" _ASM_ADC " %8,%" _ASM_PRE "1\n" \
     "\t" _ASM_ADC " %6,%" _ASM_PRE "0" \
    : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a4)), _ASM_RME ((ulong)(b4)), \
      "1" ((ulong)(a3)), _ASM_RME ((ulong)(b3)), \
      "2" ((ulong)(a2)), _ASM_RME ((ulong)(b2)), \
      "3" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "4" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define add_ssssssaaaaaaaaaaaa(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0) \
  __asm__(_ASM_ADD " %17,%" _ASM_PRE "5\n" \
     "\t" _ASM_ADC " %15,%" _ASM_PRE "4\n" \
     "\t" _ASM_ADC " %13,%" _ASM_PRE "3\n" \
     "\t" _ASM_ADC " %11,%" _ASM_PRE "2\n" \
     "\t" _ASM_ADC " %9,%" _ASM_PRE "1\n" \
     "\t" _ASM_ADC " %7,%" _ASM_PRE "0" \
    : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a5)), _ASM_RME ((ulong)(b5)), \
      "1" ((ulong)(a4)), _ASM_RME ((ulong)(b4)), \
      "2" ((ulong)(a3)), _ASM_RME ((ulong)(b3)), \
      "3" ((ulong)(a2)), _ASM_RME ((ulong)(b2)), \
      "4" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "5" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define sub_ddmmss(d1, d0, m1, m0, s1, s0) \
  __asm__(_ASM_SUB " %5,%" _ASM_PRE "1\n" \
     "\t" _ASM_SBB " %3,%" _ASM_PRE "0" \
    : "=r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "1" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

# define sub_dddmmmsss(d2, d1, d0, m2, m1, m0, s2, s1, s0) \
  __asm__(_ASM_SUB " %8,%" _ASM_PRE "2\n" \
     "\t" _ASM_SBB " %6,%" _ASM_PRE "1\n" \
     "\t" _ASM_SBB " %4,%" _ASM_PRE "0" \
    : "=r" (d2), "=&r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m2)), _ASM_RME ((ulong)(s2)), \
      "1" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "2" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

# define sub_ddddmmmmssss(d3, d2, d1, d0, m3, m2, m1, m0, s3, s2, s1, s0) \
  __asm__(_ASM_SUB " %11,%" _ASM_PRE "3\n" \
     "\t" _ASM_SBB " %9,%" _ASM_PRE "2\n" \
     "\t" _ASM_SBB " %7,%" _ASM_PRE "1\n" \
     "\t" _ASM_SBB " %5,%" _ASM_PRE "0" \
    : "=r" (d3), "=&r" (d2), "=&r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m3)), _ASM_RME ((ulong)(s3)), \
      "1" ((ulong)(m2)), _ASM_RME ((ulong)(s2)), \
      "2" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "3" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

# define sub_dddddmmmmmsssss(d4, d3, d2, d1, d0, m4, m3, m2, m1, m0, s4, s3, s2, s1, s0) \
  __asm__(_ASM_SUB " %14,%" _ASM_PRE "4\n" \
     "\t" _ASM_SBB " %12,%" _ASM_PRE "3\n" \
     "\t" _ASM_SBB " %10,%" _ASM_PRE "2\n" \
     "\t" _ASM_SBB " %8,%" _ASM_PRE "1\n" \
     "\t" _ASM_SBB " %6,%" _ASM_PRE "0" \
    : "=r" (d4), "=&r" (d3), "=&r" (d2), "=&r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m4)), _ASM_RME ((ulong)(s4)), \
      "1" ((ulong)(m3)), _ASM_RME ((ulong)(s3)), \
      "2" ((ulong)(m2)), _ASM_RME ((ulong)(s2)), \
      "3" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "4" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

# define sub_ddddddmmmmmmssssss(d5, d4, d3, d2, d1, d0, m5, m4, m3, m2, m1, m0, s5, s4, s3, s2, s1, s0) \
  __asm__(_ASM_SUB " %17,%" _ASM_PRE "5\n" \
     "\t" _ASM_SBB " %15,%" _ASM_PRE "4\n" \
     "\t" _ASM_SBB " %13,%" _ASM_PRE "3\n" \
     "\t" _ASM_SBB " %11,%" _ASM_PRE "2\n" \
     "\t" _ASM_SBB " %9,%" _ASM_PRE "1\n" \
     "\t" _ASM_SBB " %7,%" _ASM_PRE "0" \
    : "=r" (d5), "=&r" (d4), "=&r" (d3), "=&r" (d2), "=&r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m5)), _ASM_RME ((ulong)(s5)), \
      "1" ((ulong)(m4)), _ASM_RME ((ulong)(s4)), \
      "2" ((ulong)(m3)), _ASM_RME ((ulong)(s3)), \
      "3" ((ulong)(m2)), _ASM_RME ((ulong)(s2)), \
      "4" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "5" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

/* x86 does not have enough registers */
# if FLINT_BITS == 64 && defined (__amd64__)

# define add_sssssssaaaaaaaaaaaaaa(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0) \
  __asm__(_ASM_ADD " %20,%" _ASM_PRE "6\n" \
     "\t" _ASM_ADC " %18,%" _ASM_PRE "5\n" \
     "\t" _ASM_ADC " %16,%" _ASM_PRE "4\n" \
     "\t" _ASM_ADC " %14,%" _ASM_PRE "3\n" \
     "\t" _ASM_ADC " %12,%" _ASM_PRE "2\n" \
     "\t" _ASM_ADC " %10,%" _ASM_PRE "1\n" \
     "\t" _ASM_ADC " %8,%" _ASM_PRE "0" \
    : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a6)), _ASM_RME ((ulong)(b6)), \
      "1" ((ulong)(a5)), _ASM_RME ((ulong)(b5)), \
      "2" ((ulong)(a4)), _ASM_RME ((ulong)(b4)), \
      "3" ((ulong)(a3)), _ASM_RME ((ulong)(b3)), \
      "4" ((ulong)(a2)), _ASM_RME ((ulong)(b2)), \
      "5" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "6" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define add_ssssssssaaaaaaaaaaaaaaaa(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) \
  __asm__(_ASM_ADD " %23,%" _ASM_PRE "7\n" \
     "\t" _ASM_ADC " %21,%" _ASM_PRE "6\n" \
     "\t" _ASM_ADC " %19,%" _ASM_PRE "5\n" \
     "\t" _ASM_ADC " %17,%" _ASM_PRE "4\n" \
     "\t" _ASM_ADC " %15,%" _ASM_PRE "3\n" \
     "\t" _ASM_ADC " %13,%" _ASM_PRE "2\n" \
     "\t" _ASM_ADC " %11,%" _ASM_PRE "1\n" \
     "\t" _ASM_ADC " %9,%" _ASM_PRE "0" \
    : "=r" (s7), "=&r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a7)), _ASM_RME ((ulong)(b7)), \
      "1" ((ulong)(a6)), _ASM_RME ((ulong)(b6)), \
      "2" ((ulong)(a5)), _ASM_RME ((ulong)(b5)), \
      "3" ((ulong)(a4)), _ASM_RME ((ulong)(b4)), \
      "4" ((ulong)(a3)), _ASM_RME ((ulong)(b3)), \
      "5" ((ulong)(a2)), _ASM_RME ((ulong)(b2)), \
      "6" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "7" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define sub_dddddddmmmmmmmsssssss(d6, d5, d4, d3, d2, d1, d0, m6, m5, m4, m3, m2, m1, m0, s6, s5, s4, s3, s2, s1, s0) \
  __asm__(_ASM_SUB " %20,%" _ASM_PRE "6\n" \
     "\t" _ASM_SBB " %18,%" _ASM_PRE "5\n" \
     "\t" _ASM_SBB " %16,%" _ASM_PRE "4\n" \
     "\t" _ASM_SBB " %14,%" _ASM_PRE "3\n" \
     "\t" _ASM_SBB " %12,%" _ASM_PRE "2\n" \
     "\t" _ASM_SBB " %10,%" _ASM_PRE "1\n" \
     "\t" _ASM_SBB " %8,%" _ASM_PRE "0" \
    : "=r" (d6), "=&r" (d5), "=&r" (d4), "=&r" (d3), "=&r" (d2), "=&r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m6)), _ASM_RME ((ulong)(s6)), \
      "1" ((ulong)(m5)), _ASM_RME ((ulong)(s5)), \
      "2" ((ulong)(m4)), _ASM_RME ((ulong)(s4)), \
      "3" ((ulong)(m3)), _ASM_RME ((ulong)(s3)), \
      "4" ((ulong)(m2)), _ASM_RME ((ulong)(s2)), \
      "5" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "6" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

# define sub_ddddddddmmmmmmmmssssssss(d7, d6, d5, d4, d3, d2, d1, d0, m7, m6, m5, m4, m3, m2, m1, m0, s7, s6, s5, s4, s3, s2, s1, s0) \
  __asm__(_ASM_SUB " %23,%" _ASM_PRE "7\n" \
     "\t" _ASM_SBB " %21,%" _ASM_PRE "6\n" \
     "\t" _ASM_SBB " %19,%" _ASM_PRE "5\n" \
     "\t" _ASM_SBB " %17,%" _ASM_PRE "4\n" \
     "\t" _ASM_SBB " %15,%" _ASM_PRE "3\n" \
     "\t" _ASM_SBB " %13,%" _ASM_PRE "2\n" \
     "\t" _ASM_SBB " %11,%" _ASM_PRE "1\n" \
     "\t" _ASM_SBB " %9,%" _ASM_PRE "0" \
    : "=r" (d7), "=&r" (d6), "=&r" (d5), "=&r" (d4), "=&r" (d3), "=&r" (d2), "=&r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m7)), _ASM_RME ((ulong)(s7)), \
      "1" ((ulong)(m6)), _ASM_RME ((ulong)(s6)), \
      "2" ((ulong)(m5)), _ASM_RME ((ulong)(s5)), \
      "3" ((ulong)(m4)), _ASM_RME ((ulong)(s4)), \
      "4" ((ulong)(m3)), _ASM_RME ((ulong)(s3)), \
      "5" ((ulong)(m2)), _ASM_RME ((ulong)(s2)), \
      "6" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "7" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

#endif

# if defined(__BMI2__) && defined(__amd64__)
#  define umul_ppmm(w1, w0, u, v) \
  __asm__("mulx\t%3, %q0, %q1" \
    : "=r" (w0), "=r" (w1) \
    : "%d" ((ulong)(u)), "rm" ((ulong)(v)))
# else
#  define umul_ppmm(w1, w0, u, v) \
  __asm__(_ASM_MUL " %3" \
    : "=a" (w0), "=d" (w1) \
    : "%0" ((ulong)(u)), "rm" ((ulong)(v)))
#endif

# define smul_ppmm(w1, w0, u, v) \
  __asm__(_ASM_IMUL " %3" \
    : "=a" (w0), "=d" (w1) \
    : "%0" ((ulong)(u)), "rm" ((ulong)(v)))

#elif (FLINT_BITS == 64 && defined(__aarch64__)) || (FLINT_BITS == 32 && defined(__arm__))

# define add_ssaaaa(s1, s0, a1, a0, b1, b0) \
  __asm__("adds %1,%3,%5\n" \
        "\tadc %0,%2,%4" \
    : "=r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

# define sub_ddmmss(s1, s0, a1, a0, b1, b0) \
  __asm__("subs %1,%3,%5\n" \
        "\tsbc  %0,%2,%4" \
    : "=r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

# define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  __asm__("adds %2,%5,%8\n" \
        "\tadcs %1,%4,%7\n" \
        "\tadc %0,%3,%6" \
    : "=r" (s2), "=&r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

# define sub_dddmmmsss(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  __asm__("subs %2,%5,%8\n" \
        "\tsbcs %1,%4,%7\n" \
        "\tsbc  %0,%3,%6" \
    : "=r" (s2), "=&r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

/* ARM (32-bit) only has 16 registers, and has problems with inline assembly
 * with too many registers.  See issue #2131. */
# if FLINT_BITS == 64 && defined(__aarch64__)
# define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  __asm__("adds %3,%7,%11\n" \
        "\tadcs %2,%6,%10\n" \
        "\tadcs %1,%5,%9\n" \
        "\tadc  %0,%4,%8" \
    : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

# define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  __asm__ ("subs %3,%7,%11\n\tsbcs %2,%6,%10\n\tsbcs %1,%5,%9\n\tsbc %0,%4,%8"\
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")

# define add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0) \
  __asm__ ("adds %4,%9,%14\n\tadcs %3,%8,%13\n\tadcs %2,%7,%12\n\tadcs %1,%6,%11\n\tadc %0,%5,%10"\
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a4)), "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b4)), "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")

# define sub_dddddmmmmmsssss(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0) \
  __asm__ ("subs %4,%9,%14\n\tsbcs %3,%8,%13\n\tsbcs %2,%7,%12\n\tsbcs %1,%6,%11\n\tsbc %0,%5,%10"\
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a4)), "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b4)), "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")

# define add_ssssssaaaaaaaaaaaa(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0) \
  __asm__ ("adds %5,%11,%17\n\tadcs %4,%10,%16\n\tadcs %3,%9,%15\n\tadcs %2,%8,%14\n\tadcs %1,%7,%13\n\tadc %0,%6,%12"\
       : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a5)), "r" ((ulong)(a4)), "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b5)), "r" ((ulong)(b4)), "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")

# define sub_ddddddmmmmmmssssss(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0) \
  __asm__ ("subs %5,%11,%17\n\tsbcs %4,%10,%16\n\tsbcs %3,%9,%15\n\tsbcs %2,%8,%14\n\tsbcs %1,%7,%13\n\tsbc %0,%6,%12"\
       : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a5)), "r" ((ulong)(a4)), "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b5)), "r" ((ulong)(b4)), "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")

# define add_sssssssaaaaaaaaaaaaaa(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0) \
  __asm__ ("adds %6,%13,%20\n\tadcs %5,%12,%19\n\tadcs %4,%11,%18\n\tadcs %3,%10,%17\n\tadcs %2,%9,%16\n\tadcs %1,%8,%15\n\tadc %0,%7,%14"\
       : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a6)), "r" ((ulong)(a5)), "r" ((ulong)(a4)), "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b6)), "r" ((ulong)(b5)), "r" ((ulong)(b4)), "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")

# define sub_dddddddmmmmmmmsssssss(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0) \
  __asm__ ("subs %6,%13,%20\n\tsbcs %5,%12,%19\n\tsbcs %4,%11,%18\n\tsbcs %3,%10,%17\n\tsbcs %2,%9,%16\n\tsbcs %1,%8,%15\n\tsbc %0,%7,%14"\
       : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a6)), "r" ((ulong)(a5)), "r" ((ulong)(a4)), "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b6)), "r" ((ulong)(b5)), "r" ((ulong)(b4)), "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")

# define add_ssssssssaaaaaaaaaaaaaaaa(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) \
  __asm__ ("adds %7,%15,%23\n\tadcs %6,%14,%22\n\tadcs %5,%13,%21\n\tadcs %4,%12,%20\n\tadcs %3,%11,%19\n\tadcs %2,%10,%18\n\tadcs %1,%9,%17\n\tadc %0,%8,%16"\
       : "=r" (s7), "=&r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a7)), "r" ((ulong)(a6)), "r" ((ulong)(a5)), "r" ((ulong)(a4)), "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b7)), "r" ((ulong)(b6)), "r" ((ulong)(b5)), "r" ((ulong)(b4)), "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")

# define sub_ddddddddmmmmmmmmssssssss(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) \
  __asm__ ("subs %7,%15,%23\n\tsbcs %6,%14,%22\n\tsbcs %5,%13,%21\n\tsbcs %4,%12,%20\n\tsbcs %3,%11,%19\n\tsbcs %2,%10,%18\n\tsbcs %1,%9,%17\n\tsbc %0,%8,%16"\
       : "=r" (s7), "=&r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
       : "r" ((ulong)(a7)), "r" ((ulong)(a6)), "r" ((ulong)(a5)), "r" ((ulong)(a4)), "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
         "r" ((ulong)(b7)), "r" ((ulong)(b6)), "r" ((ulong)(b5)), "r" ((ulong)(b4)), "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
       : "cc")
# endif

# if defined(__arm__)
#  define umul_ppmm(xh, xl, a, b) \
  __asm__("umull %0,%1,%2,%3" \
    : "=&r" (xl), "=&r" (xh) \
    : "r" (a), "r" (b))

#  define smul_ppmm(xh, xl, a, b) \
  __asm__("smull %0,%1,%2,%3" \
    : "=&r" (xl), "=&r" (xh) \
    : "r" (a), "r" (b))
# else
#  define umul_ppmm(xh, xl, a, b) \
  __asm__("mul %0,%2,%3\n" \
      "umulh %1,%2,%3" \
    : "=&r" (xl), "=&r" (xh) \
    : "r" (a), "r" (b))

#  define smul_ppmm(xh, xl, a, b) \
  __asm__("mul %0,%2,%3\n" \
      "smulh %1,%2,%3" \
    : "=&r" (xl), "=&r" (xh) \
    : "r" (a), "r" (b))
# endif

#endif

#endif
