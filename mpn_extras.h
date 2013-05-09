/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/******************************************************************************

 Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#ifndef MPN_EXTRAS_H
#define MPN_EXTRAS_H

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define MPN_NORM(a, an)                         \
    do {                                        \
        while ((an) != 0 && (a)[(an) - 1] == 0) \
           (an)--;                              \
    } while (0)

#define MPN_SWAP(a, an, b, bn) \
    do {                       \
        mp_ptr __t;            \
        mp_size_t __tn;        \
        __t = (a);             \
        (a) = (b);             \
        (b) = __t;             \
        __tn = (an);           \
        (an) = (bn);           \
        (bn) = __tn;           \
    } while (0)

/* Not defined in gmp.h
mp_limb_t  __gmpn_modexact_1_odd(mp_srcptr src, mp_size_t size,
                                 mp_limb_t divisor);
#define mpn_modexact_1_odd __gmpn_modexact_1_odd
*/

#ifdef mpn_modexact_1_odd
#define flint_mpn_divisible_1_p(x, xsize, d) (mpn_modexact_1_odd(x, xsize, d) == 0)
#else
static __inline__ int
flint_mpn_divisible_1_p(mp_srcptr x, mp_size_t xsize, mp_limb_t d)
{
    __mpz_struct s;
    s._mp_size = xsize;
    s._mp_d = (mp_ptr) x;
    return mpz_divisible_ui_p(&s, d);
}
#endif

static __inline__
int flint_mpn_zero_p(mp_srcptr x, mp_size_t xsize)
{
    len_t i;
    for (i = 0; i < xsize; i++)
    {
        if (x[i])
            return 0;
    }
    return 1;
}

static __inline__
mp_size_t flint_mpn_divexact_1(mp_ptr x, mp_size_t xsize, mp_limb_t d)
{
    mpn_divrem_1(x, 0, x, xsize, d);
    if (x[xsize - 1] == 0UL)
        xsize -= 1;
    return xsize;
}

void flint_mpn_debug(mp_srcptr x, mp_size_t xsize);

mp_size_t flint_mpn_remove_2exp(mp_ptr x, mp_size_t xsize, mp_bitcnt_t *bits);

mp_size_t flint_mpn_remove_power_ascending(mp_ptr x, mp_size_t xsize,
                                     mp_ptr p, mp_size_t psize, ulong *exp);

int flint_mpn_factor_trial(mp_srcptr x, mp_size_t xsize, len_t start, len_t stop);

int flint_mpn_divides(mp_ptr q, mp_srcptr array1, 
         mp_size_t limbs1, mp_srcptr arrayg, mp_size_t limbsg, mp_ptr temp);

mp_size_t flint_mpn_gcd_full(mp_ptr arrayg, 
          mp_ptr array1, mp_size_t limbs1, mp_ptr array2, mp_size_t limbs2);

/*
Macros for common operations with carry management.
Also automatically order the large operand first in mpn_mul
and mpn_add.
*/

#define MPN_SET(xx, xs, yy, ys)   \
do {                              \
    xs = ys;                      \
    flint_mpn_copyi(xx, yy, ys);        \
} while (0)

#define CARRY(cy,xx,xs) if (cy) { xx[xs++] = cy; }

#define MPN_MUL_1(xx, xs, yy, ys, dd)           \
do {                                            \
    mp_limb_t __cy;                             \
    __cy = mpn_mul_1(xx, yy, ys, dd);           \
    xs = ys;                                    \
    CARRY(__cy, xx, xs);                        \
} while (0)

#define MPN_MUL(xx, xs, yy, ys, zz, zs)         \
do {                                            \
    mp_limb_t __top;                            \
    xs = ys + zs;                               \
    if (ys >= zs)                               \
        __top = mpn_mul(xx, yy, ys, zz, zs);    \
    else                                        \
        __top = mpn_mul(xx, zz, zs, yy, ys);    \
    if (!__top)                                 \
        xs--;                                   \
} while (0)

#define MPN_ADD(xx, xs, yy, ys, zz, zs)         \
do {                                            \
    mp_limb_t __cy;                             \
    if (ys >= zs)                               \
    {                                           \
        xs = ys;                                \
        __cy = mpn_add(xx, yy, ys, zz, zs);     \
    }                                           \
    else                                        \
    {                                           \
        xs = zs;                                \
        __cy = mpn_add(xx, zz, zs, yy, ys);     \
    }                                           \
    CARRY(__cy, xx, xs);                        \
} while (0)

/*
  Note: still assumes xs >= ys
*/
#define MPN_ADDMUL_1(xx, xs, yy, ys, dd)                     \
do {                                                         \
    mp_limb_t __cy;                                          \
    __cy = mpn_addmul_1(xx, yy, ys, dd);                     \
    if (__cy)                                                \
    {                                                        \
        if (xs == ys)                                        \
            xx[xs++] = __cy;                                 \
        else                                                 \
        {                                                    \
            __cy = mpn_add_1(xx+ys, xx+ys, xs-ys, __cy);     \
            CARRY(__cy, xx, xs);                             \
        }                                                    \
    }                                                        \
} while (0)



void
flint_mpn_harmonic_odd_balanced(mp_ptr t, mp_size_t * tsize,
                          mp_ptr v, mp_size_t * vsize,
                          len_t a, len_t b, len_t n, int d);

mp_limb_t flint_mpn_preinv1(mp_limb_t d, mp_limb_t d2);

mp_limb_t flint_mpn_divrem_preinv1(mp_ptr q, mp_ptr a, 
           mp_size_t m, mp_srcptr b, mp_size_t n, mp_limb_t dinv);

#define flint_mpn_divrem21_preinv(q, a_hi, a_lo, dinv) \
   do { \
      mp_limb_t __q2, __q3, __q4; \
      umul_ppmm((q), __q2, (a_hi), (dinv)); \
      umul_ppmm(__q3, __q4, (a_lo), (dinv)); \
      add_ssaaaa((q), __q2, (q), __q2, 0, __q3); \
      add_ssaaaa((q), __q2, (q), __q2, (a_hi), (a_lo)); \
   } while (0)

void flint_mpn_mulmod_preinv1(mp_ptr r, 
        mp_srcptr a, mp_srcptr b, mp_size_t n, 
        mp_srcptr d, mp_limb_t dinv, ulong norm);

void flint_mpn_preinvn(mp_ptr dinv, mp_srcptr d, mp_size_t n);

void flint_mpn_mulmod_preinvn(mp_ptr r, 
        mp_srcptr a, mp_srcptr b, mp_size_t n, 
        mp_srcptr d, mp_srcptr dinv, ulong norm);

int flint_mpn_mulmod_2expp1_basecase(mp_ptr xp, mp_srcptr yp, mp_srcptr zp, int c,
    mp_bitcnt_t b, mp_ptr tp);

static __inline__
void flint_mpn_rrandom(mp_limb_t *rp, gmp_randstate_t state, mp_size_t n)
{
  __mpz_struct str;
  str._mp_d = rp;
  str._mp_alloc = n;
  str._mp_size =n;
  mpz_rrandomb(&str,state,FLINT_BITS*n);
}

static __inline__
void flint_mpn_urandomb(mp_limb_t *rp, gmp_randstate_t state, mp_bitcnt_t n)
{
  __mpz_struct str;
  str._mp_d = rp;
  str._mp_alloc = (n + FLINT_BITS - 1)/FLINT_BITS;
  str._mp_size = (n + FLINT_BITS - 1)/FLINT_BITS;
  mpz_rrandomb(&str,state,n);
}

#ifdef __cplusplus
}
#endif

#endif
