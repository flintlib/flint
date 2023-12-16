/*
   Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
   2004, 2005 Free Software Foundation, Inc.

   Copyright 2009, 2015, 2016 William Hart
   Copyright 2011 Fredrik Johansson
   Copyright 2023 Albin Ahlbäck

   This file is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   This file is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this file; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA.
*/

/*
   N.B: This file has been adapted from code found in GMP 4.2.1.
*/

#ifndef FLINT_LONGLONG_H
#define FLINT_LONGLONG_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__GNUC__)

/* Trailing and leading zeros */
# ifndef _LONG_LONG_LIMB
#  define flint_clz __builtin_clzl
#  define flint_ctz __builtin_ctzl
# else
#  define flint_clz __builtin_clzll
#  define flint_ctz __builtin_ctzll
# endif

/* Byte swap */
# define _FLINT_CAT_(X,Y) X##Y
# define _FLINT_CAT(X,Y) _FLINT_CAT_(X,Y)
# define byte_swap(x) do { (x) = _FLINT_CAT(__builtin_bswap, GMP_LIMB_BITS)(x); } while (0)

/* Addition, subtraction and multiplication */
# if defined(__clang__)
#  include "longlong_asm_clang.h"
# else
#  include "longlong_asm_gcc.h"
# endif

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
# define NEED_CLZ_TAB
FLINT_DLL extern const unsigned char __flint_clz_tab[128];

# define flint_clz flint_clz
static inline int flint_clz(ulong x)
{
    mp_limb_t a, xr = x;
    const unsigned int bits4 = FLINT_BITS / 4;
    if (FLINT_BITS == 32)
        a = xr < ((ulong) 1 << 2 * bits4)
        ? (xr < ((ulong) 1 << bits4) ? 1 : bits4 + 1)
        : (xr < ((ulong) 1 << 3 * bits4) ? 2 * bits4 + 1 : 3 * bits4 + 1);
    else
    {
        for (__a = FLINT_BITS - 8; __a > 0; __a -= 8)
            if (((__xr >> __a) & 0xff) != 0)
                break;
        ++__a;
    }
    return FLINT_BITS + 1 - a - __flint_clz_tab[xr >> a];
}

# define flint_ctz flint_ctz
static inline int flint_ctz(ulong x)
{
    return FLINT_BITS - 1 - flint_clz(x & -x);
}
#endif

/* Byte swap */
#if !defined(byte_swap)
# if FLINT_BITS == 64
#  define byte_swap(n) \
  do { \
      /* swap adjacent bytes */ \
      n = (((n & 0xff00ff00) >> 8) | ((n & 0x00ff00ff) << 8)); \
      /* swap adjacent words */ \
      n = ((n >> 16) | (n << 16)); \
  } while (0)
# else
#  define byte_swap(n) \
  do { \
      /* swap adjacent bytes */ \
      n = (((n & 0xff00ff00ff00ff00) >> 8) | ((n & 0x00ff00ff00ff00ff) << 8)); \
      /* swap adjacent words */ \
      n = (((n & 0xffff0000ffff0000) >> 16) | ((n & 0x0000ffff0000ffff) << 16)); \
      /* swap adjacent double words */ \
      n = ((n >> 32) | (n << 32)); \
  } while (0)
# endif
#endif

/* Addition and subtraction */
#if !defined(add_ssaaaa)
# define add_ssaaaa(s1, s0, a1, a0, b1, b0) \
  do { \
    (s0) = (a0) + (b0); \
    (s1) = (a1) + (b1) + ((ulong) (s0) < (ulong) (a0)); \
  } while (0)

# define add_sssaaaaaa(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  do { \
    (s0) = (a0) + (b0); \
    (s1) = (a1) + (b1) + ((ulong) (s0) < (ulong) (a0)); \
    (s2) = (a2) + (b2) + ((ulong) (s1) < (ulong) (a1)); \
  } while (0)

# define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  do { \
    (s0) = (a0) + (b0); \
    (s1) = (a1) + (b1) + ((ulong) (s0) < (ulong) (a0)); \
    (s2) = (a2) + (b2) + ((ulong) (s1) < (ulong) (a1)); \
    (s3) = (a3) + (b3) + ((ulong) (s2) < (ulong) (a2)); \
  } while (0)

# define sub_ddmmss(s1, s0, a1, a0, b1, b0) \
  do { \
    (s0) = (a0) - (b0); \
    (s1) = (a1) - (b1) - ((ulong) (s0) > (ulong) (a0)); \
  } while (0)

# define sub_dddmmmsss(s2, s1, s0, a2, a1, a0, b2, b1, b0) \
  do { \
    (s0) = (a0) - (b0); \
    (s1) = (a1) - (b1) - ((ulong) (s0) > (ulong) (a0)); \
    (s2) = (a2) - (b2) - ((ulong) (s1) > (ulong) (a1)); \
  } while (0)
#endif

/* Multiplication */
#if !defined(umul_ppmm)
# define umul_ppmm(w1, w0, u, v) \
  do { \
    ulong __x0, __x1, __x2, __x3; \
    ulong __ul, __vl, __uh, __vh; \
    __ul = __ll_lowpart(u); \
    __uh = __ll_highpart(v); \
    __vl = __ll_lowpart(v); \
    __vh = __ll_highpart(v); \
    __x0 = __ul * __vl; \
    __x1 = __ul * __vh; \
    __x2 = __uh * __vl; \
    __x3 = __uh * __vh; \
    __x1 += __ll_highpart (__x0); /* this can't give carry */ \
    __x1 += __x2; /* but this indeed can */ \
    if (__x1 < __x2) \
      __x3 += __ll_B; \
    (w1) = __x3 + __ll_highpart (__x1); \
    (w0) = (__x1 << FLINT_BITS / 2) + __ll_lowpart(__x0); \
  } while (0)

# define smul_ppmm(w1, w0, u, v) \
  do { \
    ulong __w1; \
    umul_ppmm (__w1, w0, u, v); \
    (w1) = __w1 - (-((u) >> (FLINT_BITS - 1)) & (v)) - (-((v) >> (FLINT_BITS - 1)) & (v)); \
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
    ulong __norm = flint_clz(d); \
    if (__norm) \
    { \
      udiv_qrnnd_int((q), (r), ((n1) << __norm) + ((n0) >> (FLINT_BITS - __norm)), (n0) << __norm, (d) << __norm); \
      (r) = ((mp_limb_t) (r) >> __norm); \
    } \
    else \
       udiv_qrnnd_int((q), (r), (n1), (n0), (d)); \
  } while (0)

# define __highbit (~(ulong) 0 ^ ((~(ulong) 0) >> 1))

# define sdiv_qrnnd(q, r, n1, n0, d) \
  do { \
    mp_limb_t __n1, __n0, __d; \
    mp_limb_t __q, __r; \
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
    umul_ppmm (_xh, _xl, di, _n2 - _nmask); \
    add_ssaaaa (_xh, _xl, _xh, _xl, _n2, _nadj); \
    _q1 = ~_xh; \
    umul_ppmm (_xh, _xl, _q1, d); \
    add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl); \
    _xh -= (d); /* xh = 0 or -1 */ \
    (r) = _xl + ((d) & _xh); \
    (q) = _xh - _q1; \
  } while (0)

#ifdef __cplusplus
}
#endif

#endif
