/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CRT_HELPERS_H
#define CRT_HELPERS_H

#if defined(__GNUC__)
# if defined(__AVX2__)
#  include <x86intrin.h>
# elif defined(__ARM_NEON)
#  include <arm_neon.h>
# endif
#elif defined(_MSC_VER)
# if defined(__AVX2__)
#  include <intrin.h>
# elif defined(_M_ARM64)
#  include <arm_neon.h>
# endif
#endif

#include "flint.h"
#include "templates.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__AVX2__)

FLINT_FORCE_INLINE unsigned char _addcarry_ulong(unsigned char cf, ulong x, ulong y, ulong* s)
{
    long long unsigned int _s;
    cf = _addcarry_u64(cf, (long long unsigned int)(x),
                           (long long unsigned int)(y),
                           &_s);
    *s = (ulong)(_s);
    return cf;
}

FLINT_FORCE_INLINE unsigned char _subborrow_ulong(unsigned char cf, ulong x, ulong y, ulong* s)
{
    long long unsigned int _s;
    cf = _subborrow_u64(cf, (long long unsigned int)(x),
                            (long long unsigned int)(y),
                           &_s);
    *s = (ulong)(_s);
    return cf;
}

#else

FLINT_FORCE_INLINE unsigned char _addcarry_ulong(unsigned char cf, ulong x, ulong y, ulong* s)
{
#if 0
    ulong cf2;
    *s = __builtin_addcl(x, y, cf, &cf2);
    return cf2;
#else
    ulong hi, lo;
    add_ssaaaa(hi, lo, 0, x, 0, y);
    add_ssaaaa(hi, lo, hi, lo, 0, (ulong) cf);
    *s = lo;
    return hi;
#endif
}

FLINT_FORCE_INLINE unsigned char _subborrow_ulong(unsigned char cf, ulong x, ulong y, ulong* s)
{
#if 0
    ulong cf2;
    *s = __builtin_subcl(x, y, cf, &cf2);
    return cf2;
#else
    ulong hi, lo;
    sub_ddmmss(hi, lo, 0, x, 0, y);
    sub_ddmmss(hi, lo, hi, lo, 0, (ulong) cf);
    *s = lo;
    return hi != 0;
#endif
}

#endif


#if 1

#if defined(__GNUC__) && defined(__AVX2__)

#define add_sssssaaaaaaaaaa(s4,s3,s2,s1,s0, a4,a3,a2,a1,a0, b4,b3,b2,b1,b0)  \
  __asm__ ("addq %14,%q4\n\tadcq %12,%q3\n\tadcq %10,%q2\n\tadcq %8,%q1\n\tadcq %6,%q0"    \
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "1"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "2"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "3"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "4"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define add_ssssssaaaaaaaaaaaa(s5,s4,s3,s2,s1,s0, a5,a4,a3,a2,a1,a0, b5,b4,b3,b2,b1,b0)  \
  __asm__ ("addq %17,%q5\nadcq %15,%q4\n\tadcq %13,%q3\n\tadcq %11,%q2\n\tadcq %9,%q1\n\tadcq %7,%q0"    \
       : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "1"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "2"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "3"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "4"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "5"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define add_sssssssaaaaaaaaaaaaaa(s6,s5,s4,s3,s2,s1,s0, a6,a5,a4,a3,a2,a1,a0, b6,b5,b4,b3,b2,b1,b0)  \
  __asm__ ("addq %20,%q6\nadcq %18,%q5\nadcq %16,%q4\n\tadcq %14,%q3\n\tadcq %12,%q2\n\tadcq %10,%q1\n\tadcq %8,%q0"    \
       : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a6)), "rme" ((mp_limb_t)(b6)),                 \
         "1"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "2"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "3"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "4"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "5"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "6"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define add_ssssssssaaaaaaaaaaaaaaaa(s7,s6,s5,s4,s3,s2,s1,s0, a7,a6,a5,a4,a3,a2,a1,a0, b7,b6,b5,b4,b3,b2,b1,b0)  \
  __asm__ ("addq %23,%q7\nadcq %21,%q6\nadcq %19,%q5\n\tadcq %17,%q4\n\tadcq %15,%q3\n\tadcq %13,%q2\n\tadcq %11,%q1\n\tadcq %9,%q0"    \
       : "=r" (s7), "=&r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a7)), "rme" ((mp_limb_t)(b7)),                 \
         "1"  ((mp_limb_t)(a6)), "rme" ((mp_limb_t)(b6)),                 \
         "2"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "3"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "4"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "5"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "6"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "7"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))


#define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)  \
  __asm__ ("subq %11,%q3\n\tsbbq %9,%q2\n\tsbbq %7,%q1\n\tsbbq %5,%q0"    \
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "1"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "2"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "3"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define sub_dddddmmmmmsssss(s4,s3,s2,s1,s0, a4,a3,a2,a1,a0, b4,b3,b2,b1,b0)  \
  __asm__ ("subq %14,%q4\n\tsbbq %12,%q3\n\tsbbq %10,%q2\n\tsbbq %8,%q1\n\tsbbq %6,%q0"    \
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "1"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "2"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "3"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "4"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define sub_ddddddmmmmmmssssss(s5,s4,s3,s2,s1,s0, a5,a4,a3,a2,a1,a0, b5,b4,b3,b2,b1,b0)  \
  __asm__ ("subq %17,%q5\nsbbq %15,%q4\n\tsbbq %13,%q3\n\tsbbq %11,%q2\n\tsbbq %9,%q1\n\tsbbq %7,%q0"    \
       : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "1"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "2"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "3"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "4"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "5"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define sub_dddddddmmmmmmmsssssss(s6,s5,s4,s3,s2,s1,s0, a6,a5,a4,a3,a2,a1,a0, b6,b5,b4,b3,b2,b1,b0)  \
  __asm__ ("subq %20,%q6\nsbbq %18,%q5\nsbbq %16,%q4\n\tsbbq %14,%q3\n\tsbbq %12,%q2\n\tsbbq %10,%q1\n\tsbbq %8,%q0"    \
       : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a6)), "rme" ((mp_limb_t)(b6)),                 \
         "1"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "2"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "3"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "4"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "5"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "6"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define sub_ddddddddmmmmmmmmssssssss(s7,s6,s5,s4,s3,s2,s1,s0, a7,a6,a5,a4,a3,a2,a1,a0, b7,b6,b5,b4,b3,b2,b1,b0)  \
  __asm__ ("subq %23,%q7\nsbbq %21,%q6\nsbbq %19,%q5\n\tsbbq %17,%q4\n\tsbbq %15,%q3\n\tsbbq %13,%q2\n\tsbbq %11,%q1\n\tsbbq %9,%q0"    \
       : "=r" (s7), "=&r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a7)), "rme" ((mp_limb_t)(b7)),                 \
         "1"  ((mp_limb_t)(a6)), "rme" ((mp_limb_t)(b6)),                 \
         "2"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "3"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "4"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "5"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "6"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "7"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#elif defined(__GNUC__) && defined(__ARM_NEON)

#define add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0)      \
  __asm__ ("adds %4,%9,%14\n\tadcs %3,%8,%13\n\tadcs %2,%7,%12\n\tadcs %1,%6,%11\n\tadc %0,%5,%10"\
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a4)), "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b4)), "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0))                        \
       : "cc")

#define add_ssssssaaaaaaaaaaaa(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0)      \
  __asm__ ("adds %5,%11,%17\n\tadcs %4,%10,%16\n\tadcs %3,%9,%15\n\tadcs %2,%8,%14\n\tadcs %1,%7,%13\n\tadc %0,%6,%12"\
       : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a5)), "r" ((mp_limb_t)(a4)), "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b5)), "r" ((mp_limb_t)(b4)), "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0)) \
       : "cc")

#define add_sssssssaaaaaaaaaaaaaa(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0)      \
  __asm__ ("adds %6,%13,%20\n\tadcs %5,%12,%19\n\tadcs %4,%11,%18\n\tadcs %3,%10,%17\n\tadcs %2,%9,%16\n\tadcs %1,%8,%15\n\tadc %0,%7,%14"\
       : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a6)), "r" ((mp_limb_t)(a5)), "r" ((mp_limb_t)(a4)), "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b6)), "r" ((mp_limb_t)(b5)), "r" ((mp_limb_t)(b4)), "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0)) \
       : "cc")

#define add_ssssssssaaaaaaaaaaaaaaaa(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0)      \
  __asm__ ("adds %7,%15,%23\n\tadcs %6,%14,%22\n\tadcs %5,%13,%21\n\tadcs %4,%12,%20\n\tadcs %3,%11,%19\n\tadcs %2,%10,%18\n\tadcs %1,%9,%17\n\tadc %0,%8,%16"\
       : "=r" (s7), "=&r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a7)), "r" ((mp_limb_t)(a6)), "r" ((mp_limb_t)(a5)), "r" ((mp_limb_t)(a4)), "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b7)), "r" ((mp_limb_t)(b6)), "r" ((mp_limb_t)(b5)), "r" ((mp_limb_t)(b4)), "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0)) \
       : "cc")


#define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)      \
  __asm__ ("subs %3,%7,%11\n\tsbcs %2,%6,%10\n\tsbcs %1,%5,%9\n\tsbc %0,%4,%8"\
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0))                        \
       : "cc")

#define sub_dddddmmmmmsssss(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0)      \
  __asm__ ("subs %4,%9,%14\n\tsbcs %3,%8,%13\n\tsbcs %2,%7,%12\n\tsbcs %1,%6,%11\n\tsbc %0,%5,%10"\
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a4)), "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b4)), "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0))                        \
       : "cc")

#define sub_ddddddmmmmmmssssss(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0)      \
  __asm__ ("subs %5,%11,%17\n\tsbcs %4,%10,%16\n\tsbcs %3,%9,%15\n\tsbcs %2,%8,%14\n\tsbcs %1,%7,%13\n\tsbc %0,%6,%12"\
       : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a5)), "r" ((mp_limb_t)(a4)), "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b5)), "r" ((mp_limb_t)(b4)), "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0))                        \
       : "cc")

#define sub_dddddddmmmmmmmsssssss(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0)      \
  __asm__ ("subs %6,%13,%20\n\tsbcs %5,%12,%19\n\tsbcs %4,%11,%18\n\tsbcs %3,%10,%17\n\tsbcs %2,%9,%16\n\tsbcs %1,%8,%15\n\tsbc %0,%7,%14"\
       : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a6)), "r" ((mp_limb_t)(a5)), "r" ((mp_limb_t)(a4)), "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b6)), "r" ((mp_limb_t)(b5)), "r" ((mp_limb_t)(b4)), "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0))                        \
       : "cc")

#define sub_ddddddddmmmmmmmmssssssss(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0)      \
  __asm__ ("subs %7,%15,%23\n\tsbcs %6,%14,%22\n\tsbcs %5,%13,%21\n\tsbcs %4,%12,%20\n\tsbcs %3,%11,%19\n\tsbcs %2,%10,%18\n\tsbcs %1,%9,%17\n\tsbc %0,%8,%16"\
       : "=r" (s7), "=&r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a7)), "r" ((mp_limb_t)(a6)), "r" ((mp_limb_t)(a5)), "r" ((mp_limb_t)(a4)), "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)), "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)), \
         "r" ((mp_limb_t)(b7)), "r" ((mp_limb_t)(b6)), "r" ((mp_limb_t)(b5)), "r" ((mp_limb_t)(b4)), "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)), "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0))                        \
       : "cc")

#elif defined(_MSC_VER) && (defined(__AVX2__) || defined(_M_ARM64))
#define add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0)         \
  do {                                                                                          \
    mp_limb_t __t0 = 0;                                                                         \
    add_ssssaaaaaaaa(__t0, s2, s1, s0, (mp_limb_t) 0, a2, a1, a0, (mp_limb_t) 0, b2, b1, b0);   \
    add_ssaaaa(s4, s3, a4, a3, b4, b3);                                                         \
    add_ssaaaa(s4, s3, s4, s3, (mp_limb_t) 0, __t0);                                            \
  } while (0)

#define add_ssssssaaaaaaaaaaaa(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0)      \
  do {                                                                                                      \
    mp_limb_t __t1 = 0;                                                                                     \
    add_sssssaaaaaaaaaa(__t1, s3, s2, s1, s0, (mp_limb_t) 0, a3, a2, a1, a0, (mp_limb_t) 0, b3, b2, b1, b0);\
    add_ssaaaa(s5, s4, a5, a4, b5, b4);                                                                     \
    add_ssaaaa(s5, s4, s5, s4, (mp_limb_t) 0, __t1);                                                        \
  } while (0)

#define add_sssssssaaaaaaaaaaaaaa(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0)       \
  do {                                                                                                                      \
    mp_limb_t __t2 = 0;                                                                                                     \
    add_ssssssaaaaaaaaaaaa(__t2, s4, s3, s2, s1, s0, (mp_limb_t) 0, a4, a3, a2, a1, a0, (mp_limb_t) 0, b4, b3, b2, b1, b0); \
    add_ssaaaa(s6, s5, a6, a5, b6, b5);                                                                                     \
    add_ssaaaa(s6, s5, s6, s5, (mp_limb_t) 0, __t2);                                                                        \
  } while (0)

#define add_ssssssssaaaaaaaaaaaaaaaa(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0)        \
  do {                                                                                                                                      \
    mp_limb_t __t3 = 0;                                                                                                                     \
    add_sssssssaaaaaaaaaaaaaa(__t3, s5, s4, s3, s2, s1, s0, (mp_limb_t) 0, a5, a4, a3, a2, a1, a0, (mp_limb_t) 0, b5, b4, b3, b2, b1, b0);  \
    add_ssaaaa(s7, s6, a7, a6, b7, b6);                                                                                                     \
    add_ssaaaa(s7, s6, s7, s6, (mp_limb_t) 0, __t3);                                                                                        \
  } while (0)

#define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)        \
  do {                                                                          \
    mp_limb_t __t1, __u1;                                                       \
    sub_dddmmmsss(__t1, s1, s0, (mp_limb_t) 0, a1, a0, (mp_limb_t) 0, b1, b0);  \
    sub_ddmmss(__u1, s2, (mp_limb_t) 0, a2, (mp_limb_t) 0, b2);                 \
    sub_ddmmss(s3, s2, (a3) - (b3), s2, -__u1, -__t1);                          \
  } while (0)

#define sub_dddddmmmmmsssss(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0)         \
  do {                                                                                          \
    mp_limb_t __t2, __u2;                                                                       \
    sub_ddddmmmmssss(__t2, s2, s1, s0, (mp_limb_t) 0, a2, a1, a0, (mp_limb_t) 0, b2, b1, b0);   \
    sub_ddmmss(__u2, s3, (mp_limb_t) 0, a3, (mp_limb_t) 0, b3);                                 \
    sub_ddmmss(s4, s3, (a4) - (b4), s3, -__u2, -__t2);                                          \
  } while (0)

#define sub_ddddddmmmmmmssssss(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0)      \
  do {                                                                                                      \
    mp_limb_t __t3, __u3;                                                                                   \
    sub_dddddmmmmmsssss(__t3, s3, s2, s1, s0, (mp_limb_t) 0, a3, a2, a1, a0, (mp_limb_t) 0, b3, b2, b1, b0);\
    sub_ddmmss(__u3, s4, (mp_limb_t) 0, a4, (mp_limb_t) 0, b4);                                             \
    sub_ddmmss(s5, s4, (a5) - (b5), s4, -__u3, -__t3);                                                      \
  } while (0)

#define sub_dddddddmmmmmmmsssssss(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0)       \
  do {                                                                                                                      \
    mp_limb_t __t4, __u4;                                                                                                   \
    sub_ddddddmmmmmmssssss(__t4, s4, s3, s2, s1, s0, (mp_limb_t) 0, a4, a3, a2, a1, a0, (mp_limb_t) 0, b4, b3, b2, b1, b0); \
    sub_ddmmss(__u4, s5, (mp_limb_t) 0, a5, (mp_limb_t) 0, b5);                                                             \
    sub_ddmmss(s6, s5, (a6) - (b6), s5, -__u4, -__t4);                                                                      \
  } while (0)

#define sub_ddddddddmmmmmmmmssssssss(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0)        \
  do {                                                                                                                                      \
    mp_limb_t __t5, __u5;                                                                                                                   \
    sub_dddddddmmmmmmmsssssss(__t5, s5, s4, s3, s2, s1, s0, (mp_limb_t) 0, a5, a4, a3, a2, a1, a0, (mp_limb_t) 0, b5, b4, b3, b2, b1, b0);  \
    sub_ddmmss(__u5, s6, (mp_limb_t) 0, a6, (mp_limb_t) 0, b6);                                                                             \
    sub_ddmmss(s7, s6, (a7) - (b7), s6, -__u5, -__t5);                                                                                      \
  } while (0)

#else
# error crt_helpers.h requires AVX2 or Neon instructions
#endif

FLINT_FORCE_INLINE void multi_add_0(ulong z[], const ulong a[])
{
}

FLINT_FORCE_INLINE void multi_add_1(ulong z[], const ulong a[])
{
    z[0] += a[0];
}

FLINT_FORCE_INLINE void multi_add_2(ulong z[], const ulong a[])
{
    add_ssaaaa(z[1],z[0],
               z[1],z[0],
               a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_3(ulong z[], const ulong a[])
{
    add_sssaaaaaa(z[2],z[1],z[0],
                  z[2],z[1],z[0],
                  a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_4(ulong z[], const ulong a[])
{
    add_ssssaaaaaaaa(z[3],z[2],z[1],z[0],
                     z[3],z[2],z[1],z[0],
                     a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_5(ulong z[], const ulong a[])
{
    add_sssssaaaaaaaaaa(z[4],z[3],z[2],z[1],z[0],
                        z[4],z[3],z[2],z[1],z[0],
                        a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_6(ulong z[], const ulong a[])
{
    add_ssssssaaaaaaaaaaaa(z[5],z[4],z[3],z[2],z[1],z[0],
                           z[5],z[4],z[3],z[2],z[1],z[0],
                           a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_7(ulong z[], const ulong a[])
{
    add_sssssssaaaaaaaaaaaaaa(z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_add_8(ulong z[], const ulong a[])
{
    add_ssssssssaaaaaaaaaaaaaaaa(z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_0(ulong z[], const ulong a[])
{
}

FLINT_FORCE_INLINE void multi_sub_1(ulong z[], const ulong a[])
{
    z[0] -= a[0];
}

FLINT_FORCE_INLINE void multi_sub_2(ulong z[], const ulong a[])
{
    sub_ddmmss(z[1],z[0],
               z[1],z[0],
               a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_3(ulong z[], const ulong a[])
{
    sub_dddmmmsss(z[2],z[1],z[0],
                  z[2],z[1],z[0],
                  a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_4(ulong z[], const ulong a[])
{
    sub_ddddmmmmssss(z[3],z[2],z[1],z[0],
                     z[3],z[2],z[1],z[0],
                     a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_5(ulong z[], const ulong a[])
{
    sub_dddddmmmmmsssss(z[4],z[3],z[2],z[1],z[0],
                        z[4],z[3],z[2],z[1],z[0],
                        a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_6(ulong z[], const ulong a[])
{
    sub_ddddddmmmmmmssssss(z[5],z[4],z[3],z[2],z[1],z[0],
                           z[5],z[4],z[3],z[2],z[1],z[0],
                           a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_7(ulong z[], const ulong a[])
{
    sub_dddddddmmmmmmmsssssss(z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_sub_8(ulong z[], const ulong a[])
{
    sub_ddddddddmmmmmmmmssssssss(z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

FLINT_FORCE_INLINE void multi_rsub_0(ulong z[], const ulong a[])
{
}

FLINT_FORCE_INLINE void multi_rsub_1(ulong z[], const ulong a[])
{
    z[0] = a[0] - z[0];
}

FLINT_FORCE_INLINE void multi_rsub_2(ulong z[], const ulong a[])
{
    sub_ddmmss(z[1],z[0],
               a[1],a[0],
               z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_3(ulong z[], const ulong a[])
{
    sub_dddmmmsss(z[2],z[1],z[0],
                  a[2],a[1],a[0],
                  z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_4(ulong z[], const ulong a[])
{
    sub_ddddmmmmssss(z[3],z[2],z[1],z[0],
                     a[3],a[2],a[1],a[0],
                     z[3],z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_5(ulong z[], const ulong a[])
{
    sub_dddddmmmmmsssss(z[4],z[3],z[2],z[1],z[0],
                        a[4],a[3],a[2],a[1],a[0],
                        z[4],z[3],z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_6(ulong z[], const ulong a[])
{
    sub_ddddddmmmmmmssssss(z[5],z[4],z[3],z[2],z[1],z[0],
                           a[5],a[4],a[3],a[2],a[1],a[0],
                           z[5],z[4],z[3],z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_7(ulong z[], const ulong a[])
{
    sub_dddddddmmmmmmmsssssss(z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              a[6],a[5],a[4],a[3],a[2],a[1],a[0],
                              z[6],z[5],z[4],z[3],z[2],z[1],z[0]);
}

FLINT_FORCE_INLINE void multi_rsub_8(ulong z[], const ulong a[])
{
    sub_ddddddddmmmmmmmmssssssss(z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                                 a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0],
                                 z[7],z[6],z[5],z[4],z[3],z[2],z[1],z[0]);
}

#else


#define DEFINE_IT(n) \
FLINT_FORCE_INLINE void CAT(multi_add, n)(ulong z[], const ulong a[]) \
{ \
    unsigned char cf = 0; \
    for (ulong i = 0; i < n; i++) \
        cf = _addcarry_ulong(cf, z[i], a[i], &z[i]); \
}

DEFINE_IT(0)
DEFINE_IT(1)
DEFINE_IT(2)
DEFINE_IT(3)
DEFINE_IT(4)
DEFINE_IT(5)
DEFINE_IT(6)
DEFINE_IT(7)
DEFINE_IT(8)
#undef DEFINE_IT

#define DEFINE_IT(n) \
FLINT_FORCE_INLINE void CAT(multi_sub, n)(ulong z[], const ulong a[]) \
{ \
    unsigned char cf = 0; \
    for (ulong i = 0; i < n; i++) \
        cf = _subborrow_ulong(cf, z[i], a[i], &z[i]); \
}

DEFINE_IT(0)
DEFINE_IT(1)
DEFINE_IT(2)
DEFINE_IT(3)
DEFINE_IT(4)
DEFINE_IT(5)
DEFINE_IT(6)
DEFINE_IT(7)
#undef DEFINE_IT

#endif

#ifdef __GNUC__
FLINT_FORCE_INLINE void _mul(ulong* hi, ulong* lo, ulong y, ulong x)
{
    __uint128_t p = ((__uint128_t) x) * ((__uint128_t) y);
    *lo = (ulong) (p);
    *hi = (ulong) (p >> 64);
}

FLINT_FORCE_INLINE void _madd(ulong* hi, ulong* lo, ulong y, ulong x)
{
    __uint128_t p = ((__uint128_t) *lo) | (((__uint128_t) *hi) << 64);
    p += ((__uint128_t) x) * ((__uint128_t) y);
    *lo = (ulong) (p);
    *hi = (ulong) (p >> 64);
}
#else
FLINT_FORCE_INLINE void _mul(ulong* hi, ulong* lo, ulong y, ulong x)
{
    ulong r1, r0;
    umul_ppmm(r1, r0, x, y);
    *lo = r0;
    *hi = r1;
}

FLINT_FORCE_INLINE void _madd(ulong* hi, ulong* lo, ulong y, ulong x)
{
    ulong r1, r0;
    umul_ppmm(r1, r0, x, y);
    add_ssaaaa(*hi, *lo, r1, r0, *hi, *lo);
}
#endif

#define DEFINE_IT(n, m) \
FLINT_FORCE_INLINE void CAT3(_big_mul, n, m)(ulong r[], ulong t[], ulong C[], ulong y) \
{ \
    for (ulong k = 0; k < n; k += 2) \
    { \
        if (k + 1 < n) \
        { \
            FLINT_ASSERT(k < m); \
            _mul(&r[k+1],&r[k+0], C[k+0], y); \
        } \
        else \
        { \
            FLINT_ASSERT(k + 1 == n); \
            if (k < m) \
                r[k+0] = C[k+0]*y; \
            else \
                r[k+0] = 0; \
        } \
 \
        if (k + 2 < n) \
        { \
            FLINT_ASSERT(k + 1 < m); \
            _mul(&t[k+2],&t[k+1], C[k+1], y); \
        } \
        else if (k + 1 < n) \
        { \
            if (k + 1 < m) \
                t[k+1] = C[k+1]*y; \
            else \
                t[k+1] = 0; \
        } \
    } \
} \
FLINT_FORCE_INLINE void CAT3(_big_addmul, n, m)(ulong r[], ulong t[], ulong C[], ulong y) \
{ \
    for (ulong k = 0; k < n; k += 2) \
    { \
        if (k + 1 < n) \
        { \
            FLINT_ASSERT(k < m); \
            _madd(&r[k+1],&r[k+0], C[k+0], y); \
        } \
        else \
        { \
            FLINT_ASSERT(k + 1 == n); \
            if (k < m) \
                r[k+0] += C[k+0]*y; \
        } \
 \
        if (k + 2 < n) \
        { \
            FLINT_ASSERT(k + 1 < m); \
            _madd(&t[k+2],&t[k+1], C[k+1], y); \
        } \
        else if (k + 1 < n) \
        { \
            if (k + 1 < m) \
                t[k+1] += C[k+1]*y; \
        } \
    } \
}

DEFINE_IT(1, 0)
DEFINE_IT(2, 1)
DEFINE_IT(3, 2)
DEFINE_IT(4, 3)
DEFINE_IT(4, 4)
DEFINE_IT(5, 4)
DEFINE_IT(6, 5)
DEFINE_IT(7, 6)
#undef DEFINE_IT



#define DEFINE_IT(n, n_minus_1) \
FLINT_FORCE_INLINE void CAT(_reduce_big_sum, n)(ulong r[], ulong t[], const ulong* limit) \
{ \
    CAT(multi_add, n_minus_1)(r+1, t+1); \
check: \
    for (ulong k = n; k > 1; k--) \
    { \
        if (FLINT_LIKELY(r[k-1] > limit[k-1])) \
            goto sub; \
        if (r[k-1] < limit[k-1]) \
            return; \
    } \
    if (r[0] < limit[0]) \
        return; \
sub: \
    CAT(multi_sub, n)(r, limit); \
    goto check; \
}

DEFINE_IT(1, 0)
DEFINE_IT(2, 1)
DEFINE_IT(3, 2)
DEFINE_IT(4, 3)
DEFINE_IT(5, 4)
DEFINE_IT(6, 5)
DEFINE_IT(7, 6)
#undef DEFINE_IT

#ifdef __cplusplus
}
#endif

#endif /* CRT_HELPERS_H */
