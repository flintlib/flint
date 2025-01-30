/*
    Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
    2004, 2005 Free Software Foundation, Inc.

    Copyright (C) 2009, 2015, 2016 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    Contains code from GNU MP Library.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_LONGLONG_H
#define FLINT_LONGLONG_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

/* FIXME: Bug in Clang 16+ */
#if __clang_major__ >= 16
# define CLANG_volatile volatile
#else
# define CLANG_volatile
#endif

#if defined(__GNUC__)

/* Trailing and leading zeros */
# if FLINT_LONG_LONG
#  define flint_clz __builtin_clzll
#  define flint_ctz __builtin_ctzll
# else
#  define flint_clz __builtin_clzl
#  define flint_ctz __builtin_ctzl
# endif

/* Addition, subtraction and multiplication */
# if defined(__clang__)
#  include "longlong_asm_clang.h"
# else
#  include "longlong_asm_gcc.h"
# endif
# include "longlong_asm_gnu.h"

/* Division */
# include "longlong_div_gnu.h"

#elif defined(_MSC_VER)

# if defined(_M_X64) || defined(_M_IX86)
#  include "longlong_msc_x86.h"
# elif defined(_M_ARM64)
#  include "longlong_msc_arm64.h"
# endif

#endif

/* Generics ******************************************************************/

#define __ll_B ((ulong) 1 << (FLINT_BITS / 2))
#define __ll_lowpart(t) ((ulong) (t) & (__ll_B - 1))
#define __ll_highpart(t) ((ulong) (t) >> (FLINT_BITS / 2))

/* Trailing and leading zeros */
#if !defined(flint_ctz)
# define NEED_CLZ_TAB 1
FLINT_DLL extern const unsigned char __flint_clz_tab[128];

# define flint_clz flint_clz
static inline int flint_clz(ulong x)
{
    ulong a, xr = x;
    const unsigned int bits4 = FLINT_BITS / 4;
    if (FLINT_BITS == 32)
        a = xr < ((ulong) 1 << 2 * bits4)
        ? (xr < ((ulong) 1 << bits4) ? 1 : bits4 + 1)
        : (xr < ((ulong) 1 << 3 * bits4) ? 2 * bits4 + 1 : 3 * bits4 + 1);
    else
    {
        for (a = FLINT_BITS - 8; a > 0; a -= 8)
            if (((xr >> a) & 0xff) != 0)
                break;
        ++a;
    }
    return FLINT_BITS + 1 - a - __flint_clz_tab[xr >> a];
}

# define flint_ctz flint_ctz
static inline int flint_ctz(ulong x)
{
    return FLINT_BITS - 1 - flint_clz(x & -x);
}
#endif

/* Beware when using the unsigned return value in signed arithmetic */
FLINT_FORCE_INLINE
flint_bitcnt_t FLINT_BIT_COUNT(ulong x)
{
    flint_bitcnt_t zeros = FLINT_BITS;
    if (x) zeros = flint_clz(x);
    return FLINT_BITS - zeros;
}

#define FLINT_FLOG2(k) (FLINT_BIT_COUNT(k) - 1)
#define FLINT_CLOG2(k) FLINT_BIT_COUNT((k) - 1)

/* Addition and subtraction */
#if !defined(add_ssaaaa)
# define add_ssaaaa(s1, s0, a1, a0, b1, b0) \
  do { \
    ulong __t0 = (a0); \
    (s0) = (a0) + (b0); \
    (s1) = (a1) + (b1) + ((ulong) (s0) < __t0); \
  } while (0)

# define sub_ddmmss(s1, s0, a1, a0, b1, b0) \
  do { \
    ulong __t0 = (a0); \
    (s0) = (a0) - (b0); \
    (s1) = (a1) - (b1) - ((ulong) (s0) > __t0); \
  } while (0)

# define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  do { \
    ulong __t1, __t2; \
    add_ssaaaa(__t1, s0, (ulong) 0, a0, (ulong) 0, b0); \
    add_ssaaaa(__t2, s1, (ulong) 0, a1, (ulong) 0, b1); \
    add_ssaaaa(s2, s1, (a2) + (b2), s1, __t2, __t1); \
  } while (0)

# define sub_dddmmmsss(d2, d1, d0, m2, m1, m0, s2, s1, s0) \
  do { \
    ulong __t1, __t2; \
    sub_ddmmss(__t1, d0, (ulong) 0, m0, (ulong) 0, s0); \
    sub_ddmmss(__t2, d1, (ulong) 0, m1, (ulong) 0, s1); \
    sub_ddmmss(d2, d1, (m2) - (s2), d1, -__t2, -__t1); \
  } while (0)
#endif

#if !defined(add_ssssaaaaaaaa)
# define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  do { \
    ulong __t3; \
    add_sssaaaaaa(__t3, s1, s0, (ulong) 0, a1, a0, (ulong) 0, b1, b0); \
    add_ssaaaa(s3, s2, a3, a2, b3, b2); \
    add_ssaaaa(s3, s2, s3, s2, (ulong) 0, __t3); \
  } while (0)

# define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  do { \
    ulong __t3, __t4; \
    sub_dddmmmsss(__t3, s1, s0, (ulong) 0, a1, a0, (ulong) 0, b1, b0); \
    sub_ddmmss(__t4, s2, (ulong) 0, a2, (ulong) 0, b2); \
    sub_ddmmss(s3, s2, (a3) - (b3), s2, -__t4, -__t3); \
  } while (0)

# define add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0) \
  do { \
    ulong __t4 = 0; \
    add_ssssaaaaaaaa(__t4, s2, s1, s0, (ulong) 0, a2, a1, a0, (ulong) 0, b2, b1, b0); \
    add_ssaaaa(s4, s3, a4, a3, b4, b3); \
    add_ssaaaa(s4, s3, s4, s3, (ulong) 0, __t4); \
  } while (0)

# define sub_dddddmmmmmsssss(s4, s3, s2, s1, s0, a4, a3, a2, a1, a0, b4, b3, b2, b1, b0) \
  do { \
    ulong __t5, __t6; \
    sub_ddddmmmmssss(__t5, s2, s1, s0, (ulong) 0, a2, a1, a0, (ulong) 0, b2, b1, b0); \
    sub_ddmmss(__t6, s3, (ulong) 0, a3, (ulong) 0, b3); \
    sub_ddmmss(s4, s3, (a4) - (b4), s3, -__t6, -__t5); \
  } while (0)

# define add_ssssssaaaaaaaaaaaa(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0) \
  do { \
    ulong __t5 = 0; \
    add_sssssaaaaaaaaaa(__t5, s3, s2, s1, s0, (ulong) 0, a3, a2, a1, a0, (ulong) 0, b3, b2, b1, b0); \
    add_ssaaaa(s5, s4, a5, a4, b5, b4); \
    add_ssaaaa(s5, s4, s5, s4, (ulong) 0, __t5); \
  } while (0)

# define sub_ddddddmmmmmmssssss(s5, s4, s3, s2, s1, s0, a5, a4, a3, a2, a1, a0, b5, b4, b3, b2, b1, b0) \
  do { \
    ulong __t7, __t8; \
    sub_dddddmmmmmsssss(__t7, s3, s2, s1, s0, (ulong) 0, a3, a2, a1, a0, (ulong) 0, b3, b2, b1, b0);\
    sub_ddmmss(__t8, s4, (ulong) 0, a4, (ulong) 0, b4); \
    sub_ddmmss(s5, s4, (a5) - (b5), s4, -__t8, -__t7); \
  } while (0)
#endif

/* extra wide variants might not have inline asm if there are not enough registers */
#if !defined(add_sssssssaaaaaaaaaaaaaa)
# define add_sssssssaaaaaaaaaaaaaa(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0) \
  do { \
    ulong __t6 = 0; \
    add_ssssssaaaaaaaaaaaa(__t6, s4, s3, s2, s1, s0, (ulong) 0, a4, a3, a2, a1, a0, (ulong) 0, b4, b3, b2, b1, b0); \
    add_ssaaaa(s6, s5, a6, a5, b6, b5); \
    add_ssaaaa(s6, s5, s6, s5, (ulong) 0, __t6); \
  } while (0)

# define add_ssssssssaaaaaaaaaaaaaaaa(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) \
  do { \
    ulong __t7 = 0; \
    add_sssssssaaaaaaaaaaaaaa(__t7, s5, s4, s3, s2, s1, s0, (ulong) 0, a5, a4, a3, a2, a1, a0, (ulong) 0, b5, b4, b3, b2, b1, b0); \
    add_ssaaaa(s7, s6, a7, a6, b7, b6); \
    add_ssaaaa(s7, s6, s7, s6, (ulong) 0, __t7); \
  } while (0)

# define sub_dddddddmmmmmmmsssssss(s6, s5, s4, s3, s2, s1, s0, a6, a5, a4, a3, a2, a1, a0, b6, b5, b4, b3, b2, b1, b0) \
  do { \
    ulong __t9, __t10; \
    sub_ddddddmmmmmmssssss(__t9, s4, s3, s2, s1, s0, (ulong) 0, a4, a3, a2, a1, a0, (ulong) 0, b4, b3, b2, b1, b0); \
    sub_ddmmss(__t10, s5, (ulong) 0, a5, (ulong) 0, b5); \
    sub_ddmmss(s6, s5, (a6) - (b6), s5, -__t10, -__t9); \
  } while (0)

# define sub_ddddddddmmmmmmmmssssssss(s7, s6, s5, s4, s3, s2, s1, s0, a7, a6, a5, a4, a3, a2, a1, a0, b7, b6, b5, b4, b3, b2, b1, b0) \
  do { \
    ulong __t11, __t12; \
    sub_dddddddmmmmmmmsssssss(__t11, s5, s4, s3, s2, s1, s0, (ulong) 0, a5, a4, a3, a2, a1, a0, (ulong) 0, b5, b4, b3, b2, b1, b0); \
    sub_ddmmss(__t12, s6, (ulong) 0, a6, (ulong) 0, b6); \
    sub_ddmmss(s7, s6, (a7) - (b7), s6, -__t12, -__t11); \
  } while (0)
#endif


#if !defined(MPN_INCR_U)
# if FLINT_WANT_ASSERT
#  define MPN_INCR_U(ptr, size, incr) \
  do { \
    ulong __cy = mpn_add_1(ptr, ptr, size, incr); \
    FLINT_ASSERT(__cy == 0); \
  } while (0)
#  define MPN_DECR_U(ptr, size, incr) \
  do { \
    ulong __cy = mpn_sub_1(ptr, ptr, size, incr); \
    FLINT_ASSERT(__cy == 0); \
  } while (0)
# else
#  define MPN_INCR_U(ptr, size, incr) mpn_add_1(ptr, ptr, size, incr)
#  define MPN_DECR_U(ptr, size, incr) mpn_sub_1(ptr, ptr, size, incr)
# endif
#endif

/* Multiplication */
#if !defined(umul_ppmm)
# define umul_ppmm(w1, w0, u, v) \
  do { \
    ulong __x0, __x1, __x2, __x3; \
    ulong __ul, __vl, __uh, __vh; \
    __ul = __ll_lowpart(u); \
    __uh = __ll_highpart(u); \
    __vl = __ll_lowpart(v); \
    __vh = __ll_highpart(v); \
    __x0 = __ul * __vl; \
    __x1 = __ul * __vh; \
    __x2 = __uh * __vl; \
    __x3 = __uh * __vh; \
    __x1 += __ll_highpart(__x0); /* this can't give carry */ \
    __x1 += __x2; /* but this indeed can */ \
    if (__x1 < __x2) \
      __x3 += __ll_B; \
    (w1) = __x3 + __ll_highpart(__x1); \
    (w0) = (__x1 << (FLINT_BITS / 2)) + __ll_lowpart(__x0); \
  } while (0)

# define smul_ppmm(w1, w0, u, v) \
  do { \
    ulong __w1, __u = (u), __v = (v); \
    umul_ppmm(__w1, w0, u, v); \
    (w1) = __w1 - (-(__u >> (FLINT_BITS - 1)) & __v) - (-(__v >> (FLINT_BITS - 1)) & __u); \
  } while (0)
#endif

/* Division */
#if !defined(udiv_qrnnd)
# define udiv_qrnnd_int(q, r, n1, n0, d) \
  do { \
    ulong __d1, __d0, __q1, __q0, __r1, __r0, __m; \
    FLINT_ASSERT((d) != 0); \
    FLINT_ASSERT((n1) < (d)); \
    __d1 = __ll_highpart (d); \
    __d0 = __ll_lowpart (d); \
    __q1 = (n1) / __d1; \
    __r1 = (n1) - __q1 * __d1; \
    __m = __q1 * __d0; \
    __r1 = __r1 * __ll_B | __ll_highpart (n0); \
    if (__r1 < __m) \
    { \
      __q1--, __r1 += (d); \
      if (__r1 >= (d)) /* i.e. we didn't get carry when adding to __r1 */ \
        if (__r1 < __m) \
          __q1--, __r1 += (d); \
    } \
    __r1 -= __m; \
    __q0 = __r1 / __d1; \
    __r0 = __r1  - __q0 * __d1; \
    __m = __q0 * __d0; \
    __r0 = __r0 * __ll_B | __ll_lowpart (n0); \
    if (__r0 < __m) \
    { \
      __q0--, __r0 += (d); \
      if (__r0 >= (d)) \
        if (__r0 < __m) \
          __q0--, __r0 += (d); \
    } \
    __r0 -= __m; \
    (q) = __q1 * __ll_B | __q0; \
    (r) = __r0; \
  } while (0)

# define udiv_qrnnd(q, r, n1, n0, d) \
  do { \
    CLANG_volatile int __norm = flint_clz(d); \
    if (__norm != 0) \
    { \
      udiv_qrnnd_int((q), (r), ((n1) << __norm) + ((n0) >> (FLINT_BITS - __norm)), (n0) << __norm, (d) << __norm); \
      (r) = ((ulong) (r) >> __norm); \
    } \
    else \
       udiv_qrnnd_int((q), (r), (n1), (n0), (d)); \
  } while (0)

# define __highbit (~(ulong) 0 ^ ((~(ulong) 0) >> 1))

# define sdiv_qrnnd(q, r, n1, n0, d) \
  do { \
    ulong __n1, __n0, __d; \
    ulong __q, __r; \
    unsigned int __sgn_n = 0, __sgn_d = 0; \
    if ((n1) & __highbit) \
    { \
      __n0 = -(n0); \
      __n1 = ~(n1) + (__n0 == 0); \
      __sgn_n = ~__sgn_n; \
    } \
    else \
    { \
      __n0 = (n0); \
      __n1 = (n1); \
    } \
    if ((d) & __highbit) \
    { \
        __d = -(d); \
        __sgn_d = ~__sgn_d; \
    } \
    else \
    { \
        __d = (d); \
    } \
    udiv_qrnnd(__q, __r, __n1, __n0, __d); \
    q = (__sgn_n == __sgn_d) ? __q : -__q; \
    r = (__sgn_n == 0) ? __r : -__r; \
  } while (0)
#endif

#define udiv_qrnnd_preinv(q, r, nh, nl, d, di) \
  do { \
    ulong _n2, _n10, _nmask, _nadj, _q1; \
    ulong _xh, _xl; \
    _n2 = (nh); \
    _n10 = (nl); \
    _nmask = (slong) (_n10) >> (FLINT_BITS - 1); \
    _nadj = _n10 + (_nmask & (d)); \
    umul_ppmm(_xh, _xl, di, _n2 - _nmask); \
    add_ssaaaa(_xh, _xl, _xh, _xl, _n2, _nadj); \
    _q1 = ~_xh; \
    umul_ppmm(_xh, _xl, _q1, d); \
    add_ssaaaa(_xh, _xl, _xh, _xl, nh, nl); \
    _xh -= (d); /* xh = 0 or -1 */ \
    (r) = _xl + ((d) & _xh); \
    (q) = _xh - _q1; \
  } while (0)

#ifdef __cplusplus
}
#endif

#endif
