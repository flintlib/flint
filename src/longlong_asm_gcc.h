/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LONGLONG_ASM_H
#define LONGLONG_ASM_H

/* Machine specific operations */
#if defined (__amd64__) || (GMP_LIMB_BITS == 32 && (defined (__i386__) || defined (__i486__)))

# if GMP_LIMB_BITS == 64 && defined (__amd64__)
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
  __asm__(_ASM_ADD " %5,%" _ASM_PRE "1\n"   \
     "\t" _ASM_ADC " %3,%" _ASM_PRE "0"     \
    : "=r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "%1" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  __asm__(_ASM_ADD " %8,%" _ASM_PRE "2\n"   \
     "\t" _ASM_ADC " %6,%" _ASM_PRE "1\n"   \
     "\t" _ASM_ADC " %4,%" _ASM_PRE "0"     \
    : "=r" (s2), "=&r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a2)), _ASM_RME ((ulong)(b2)), \
      "1" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "2" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  __asm__(_ASM_ADD " %11,%" _ASM_PRE"3\n"   \
     "\t" _ASM_ADC " %9,%" _ASM_PRE "2\n"   \
     "\t" _ASM_ADC " %7,%" _ASM_PRE "1\n"   \
     "\t" _ASM_ADC " %5,%" _ASM_PRE "0"     \
    : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
    : "0" ((ulong)(a3)), _ASM_RME ((ulong)(b3)), \
      "1" ((ulong)(a2)), _ASM_RME ((ulong)(b2)), \
      "2" ((ulong)(a1)), _ASM_RME ((ulong)(b1)), \
      "3" ((ulong)(a0)), _ASM_RME ((ulong)(b0)))

# define sub_ddmmss(d1, d0, m1, m0, s1, s0) \
  __asm__(_ASM_SUB " %5,%" _ASM_PRE "1\n"   \
     "\t" _ASM_SBB " %3,%" _ASM_PRE "0"     \
    : "=r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "1" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

# define sub_dddmmmsss(d2, d1, d0, m2, m1, m0, s2, s1, s0) \
  __asm__(_ASM_SUB " %8,%" _ASM_PRE "2\n"   \
     "\t" _ASM_SBB " %6,%" _ASM_PRE "1\n"   \
     "\t" _ASM_SBB " %4,%" _ASM_PRE "0"     \
    : "=r" (d2), "=&r" (d1), "=&r" (d0) \
    : "0" ((ulong)(m2)), _ASM_RME ((ulong)(s2)), \
      "1" ((ulong)(m1)), _ASM_RME ((ulong)(s1)), \
      "2" ((ulong)(m0)), _ASM_RME ((ulong)(s0)))

# define umul_ppmm(w1, w0, u, v) \
  __asm__(_ASM_MUL " %3" \
    : "=a" (w0), "=d" (w1) \
    : "%0" ((ulong)(u)), "rm" ((ulong)(v)))

# define smul_ppmm(w1, w0, u, v) \
  __asm__(_ASM_IMUL " %3" \
    : "=a" (w0), "=d" (w1) \
    : "%0" ((ulong)(u)), "rm" ((ulong)(v)))

#elif (GMP_LIMB_BITS == 64 && defined(__aarch64__)) || (GMP_LIMB_BITS == 32 && defined(__arm__))

# define add_ssaaaa(s1, s0, a1, a0, b1, b0) \
  __asm__("adds %1,%3,%5\n" \
        "\tadc %0,%2,%4"    \
    : "=r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

# define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  __asm__("adds %2,%5,%8\n" \
        "\tadcs %1,%4,%7\n" \
        "\tadc %0,%3,%6"    \
    : "=r" (s2), "=&r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

# define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  __asm__("adds %3,%7,%11\n"    \
        "\tadcs %2,%6,%10\n"    \
        "\tadcs %1,%5,%9\n"     \
        "\tadc  %0,%4,%8"       \
    : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a3)), "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b3)), "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

# define sub_ddmmss(s1, s0, a1, a0, b1, b0) \
  __asm__("subs %1,%3,%5\n" \
        "\tsbc  %0,%2,%4"   \
    : "=r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

# define sub_dddmmmsss(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  __asm__("subs %2,%5,%8\n" \
        "\tsbcs %1,%4,%7\n" \
        "\tsbc  %0,%3,%6"   \
    : "=r" (s2), "=&r" (s1), "=&r" (s0) \
    : "r" ((ulong)(a2)), "r" ((ulong)(a1)), "r" ((ulong)(a0)), \
      "r" ((ulong)(b2)), "r" ((ulong)(b1)), "rI" ((ulong)(b0)) \
    : "cc")

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
