/*
    Copyright (C) 2010 Fredrik Johansson

    2x2 mul code taken from MPFR 2.3.0
    (Copyright (C) 1991-2007 Free Software Foundation, Inc.)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_EXTRAS_H
#define MPN_EXTRAS_H

#ifdef MPN_EXTRAS_INLINES_C
#define MPN_EXTRAS_INLINE
#else
#define MPN_EXTRAS_INLINE static inline
#endif

#include "flint.h"

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

#define BITS_TO_LIMBS(b) (((b) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

double
flint_mpn_get_d(mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp);

#define FLINT_MPN_MUL_2X1(r2, r1, r0, a1, a0, b0)           \
    do {                                                    \
        mp_limb_t t1;                                       \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
    } while (0)

/* todo: use optimal code here */
#define FLINT_MPN_MUL_2X2(r3, r2, r1, r0, a1, a0, b1, b0)   \
    do {                                                    \
        mp_limb_t t1, t2, t3;                               \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
        umul_ppmm(t1, t2, a0, b1);                          \
        umul_ppmm(r3, t3, a1, b1);                          \
        add_ssaaaa(r3, t1, r3, t1, 0, t3);                  \
        add_ssaaaa(r2, r1, r2, r1, t1, t2);                 \
        r3 += r2 < t1;                                      \
    } while (0)

#define FLINT_MPN_MUL_THRESHOLD 400

mp_limb_t flint_mpn_mul_large(mp_ptr r1, mp_srcptr i1, mp_size_t n1, mp_srcptr i2, mp_size_t n2);

void flint_mpn_mul_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n);

mp_limb_t flint_mpn_mul(mp_ptr z, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn);

MPN_EXTRAS_INLINE void
flint_mpn_sqr(mp_ptr z, mp_srcptr x, mp_size_t n)
{
    if (n < FLINT_MPN_MUL_THRESHOLD)
        mpn_sqr(z, x, n);
    else
        flint_mpn_mul_large(z, x, n, x, n);
}

#define FLINT_MPN_MUL_WITH_SPECIAL_CASES(_z, _x, _xn, _y, _yn) \
    if ((_xn) == (_yn)) \
    { \
        if ((_xn) == 1) \
        { \
            umul_ppmm((_z)[1], (_z)[0], (_x)[0], (_y)[0]); \
        } \
        else if ((_xn) == 2) \
        { \
            mp_limb_t __tt_x1, __tt_x0, __tt_y1, __tt_y0; \
            __tt_x0 = (_x)[0]; \
            __tt_x1 = (_x)[1]; \
            __tt_y0 = (_y)[0]; \
            __tt_y1 = (_y)[1]; \
            FLINT_MPN_MUL_2X2((_z)[3], (_z)[2], (_z)[1], (_z)[0], __tt_x1, __tt_x0, __tt_y1, __tt_y0); \
        } \
        else if ((_x) == (_y)) \
            flint_mpn_sqr((_z), (_x), (_xn)); \
        else \
            flint_mpn_mul_n((_z), (_x), (_y), (_xn)); \
    } \
    else if ((_xn) > (_yn)) \
    { \
        if ((_yn) == 1) \
            (_z)[(_xn) + (_yn) - 1] = mpn_mul_1((_z), (_x), (_xn), (_y)[0]); \
        else \
            flint_mpn_mul((_z), (_x), (_xn), (_y), (_yn)); \
    } \
    else \
    { \
        if ((_xn) == 1) \
            (_z)[(_xn) + (_yn) - 1] = mpn_mul_1((_z), (_y), (_yn), (_x)[0]); \
        else \
            flint_mpn_mul((_z), (_y), (_yn), (_x), (_xn)); \
    }

/*
    return the high limb of a two limb left shift by n < GMP_LIMB_BITS bits.
    Note: if GMP_NAIL_BITS != 0, the rest of flint is already broken anyways.
*/
#define MPN_LEFT_SHIFT_HI(hi, lo, n)                                \
    ((n) > 0 ? (((hi) << (n)) | ((lo) >> (GMP_LIMB_BITS - (n))))    \
             : (hi))

#define MPN_RIGHT_SHIFT_LOW(hi, lo, n)                                \
    ((n) > 0 ? (((lo) >> (n)) | ((hi) << (GMP_LIMB_BITS - (n))))    \
             : (lo))

#ifdef FLINT_HAVE_MPN_MODEXACT_1_ODD

# define mpn_modexact_1_odd __gmpn_modexact_1_odd
mp_limb_t mpn_modexact_1_odd(mp_srcptr, mp_size_t, mp_limb_t);

MPN_EXTRAS_INLINE int
flint_mpn_divisible_1_odd(mp_srcptr x, mp_size_t xsize, mp_limb_t d)
{
    return mpn_modexact_1_odd(x, xsize, d) == 0;
}

#else

# include "gmpcompat.h"

MPN_EXTRAS_INLINE int
flint_mpn_divisible_1_odd(mp_srcptr x, mp_size_t xsize, mp_limb_t d)
{
    __mpz_struct s;
    s._mp_size = xsize;
    s._mp_d = (mp_ptr) x;
    return flint_mpz_divisible_ui_p(&s, d);
}

#endif

static inline void
mpn_tdiv_q(mp_ptr qp, mp_srcptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn)
{
    mp_ptr _scratch;
    TMP_INIT;
    TMP_START;
    _scratch = (mp_ptr) TMP_ALLOC(dn * sizeof(mp_limb_t));
    mpn_tdiv_qr(qp, _scratch, 0, np, nn, dp, dn);
    TMP_END;
}

MPN_EXTRAS_INLINE
int flint_mpn_zero_p(mp_srcptr x, mp_size_t xsize)
{
    slong i;
    for (i = 0; i < xsize; i++)
    {
        if (x[i])
            return 0;
    }
    return 1;
}

MPN_EXTRAS_INLINE
mp_size_t flint_mpn_divexact_1(mp_ptr x, mp_size_t xsize, mp_limb_t d)
{
    mpn_divrem_1(x, 0, x, xsize, d);
    if (x[xsize - 1] == UWORD(0))
        xsize -= 1;
    return xsize;
}

void flint_mpn_debug(mp_srcptr x, mp_size_t xsize);

mp_size_t flint_mpn_remove_2exp(mp_ptr x, mp_size_t xsize,
		                                      flint_bitcnt_t *bits);

mp_size_t flint_mpn_remove_power_ascending(mp_ptr x,
		    mp_size_t xsize, mp_ptr p, mp_size_t psize, ulong *exp);

int flint_mpn_factor_trial(mp_srcptr x, mp_size_t xsize,
		                                   slong start, slong stop);

int flint_mpn_factor_trial_tree(slong * factors,
                            mp_srcptr x, mp_size_t xsize, slong num_primes);

mp_size_t flint_mpn_fmms1(mp_ptr y, mp_limb_t a1, mp_srcptr x1,
                                      mp_limb_t a2, mp_srcptr x2, mp_size_t n);

int flint_mpn_divides(mp_ptr q, mp_srcptr array1,
         mp_size_t limbs1, mp_srcptr arrayg, mp_size_t limbsg, mp_ptr temp);

mp_size_t flint_mpn_gcd_full2(mp_ptr arrayg,
		                 mp_srcptr array1, mp_size_t limbs1,
			   mp_srcptr array2, mp_size_t limbs2, mp_ptr temp);

mp_size_t flint_mpn_gcd_full(mp_ptr arrayg,
    mp_srcptr array1, mp_size_t limbs1, mp_srcptr array2, mp_size_t limbs2);

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

void flint_mpn_mod_preinvn(mp_ptr r, mp_srcptr a, mp_size_t m,
                                     mp_srcptr d, mp_size_t n, mp_srcptr dinv);

mp_limb_t flint_mpn_divrem_preinvn(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t m,
                                     mp_srcptr d, mp_size_t n, mp_srcptr dinv);

void flint_mpn_mulmod_preinvn(mp_ptr r,
        mp_srcptr a, mp_srcptr b, mp_size_t n,
        mp_srcptr d, mp_srcptr dinv, ulong norm);

int flint_mpn_mulmod_2expp1_basecase(mp_ptr xp, mp_srcptr yp, mp_srcptr zp,
    int c, flint_bitcnt_t b, mp_ptr tp);

MPN_EXTRAS_INLINE
void flint_mpn_rrandom(mp_limb_t *rp, gmp_randstate_t state, mp_size_t n)
{
  __mpz_struct str;
  str._mp_d = rp;
  str._mp_alloc = n;
  str._mp_size =n;
  mpz_rrandomb(&str,state,FLINT_BITS*n);
}

MPN_EXTRAS_INLINE
void flint_mpn_urandomb(mp_limb_t *rp, gmp_randstate_t state, flint_bitcnt_t n)
{
  __mpz_struct str;
  str._mp_d = rp;
  str._mp_alloc = (n + FLINT_BITS - 1)/FLINT_BITS;
  str._mp_size = (n + FLINT_BITS - 1)/FLINT_BITS;
  mpz_rrandomb(&str,state,n);
}


/******************************************************************************
    Divisions where the quotient is expected to be small. All function do:
        input: n > d > 0
        output: q = n/d, r = n%d
    for various small sizes of n and d.
    Not in a function because compiler refuses to inline eudiv_qrrnndd.
    Each macro takes a prefix t for its local vars.
******************************************************************************/

#define eudiv_qrnd(q, r, n, d, t)           \
do {                                        \
    mp_limb_t t##q, t##a = n, t##b = d;     \
                                            \
    FLINT_ASSERT(t##a > t##b);              \
    FLINT_ASSERT(t##b > 0);                 \
                                            \
    t##a -= t##b;                           \
    for (t##q = 1; t##q < 5; t##q++)        \
    {                                       \
        if (t##a < t##b)                    \
            goto t##quotient_found;         \
        t##a -= t##b;                       \
    }                                       \
    t##q += t##a / t##b;                    \
    t##a = t##a % t##b;                     \
                                            \
t##quotient_found:                          \
                                            \
    q = t##q;                               \
    r = t##a;                               \
} while (0)


#define eudiv_qqrnnd(q1, q0, r0, n1, n0, d0, t)         \
do {                                                    \
    mp_limb_t t##a1 = n1, t##a0 = n0, t##b0 = d0;       \
    mp_limb_t t##q1, t##q0, t##r0, t##u;                \
                                                        \
    FLINT_ASSERT(t##a1 > 0 || t##a0 > t##b0);           \
                                                        \
    udiv_qrnnd(t##q1, t##u, 0, t##a1, t##b0);           \
    udiv_qrnnd(t##q0, t##r0, t##u, t##a0, t##b0);       \
                                                        \
    q1 = t##q1;                                         \
    q0 = t##q0;                                         \
    r0 = t##r0;                                         \
} while (0)

/* d must be normalized, i.e. d1 != 0 */
#define eudiv_qrrnndd(q0, r1, r0, n1, n0, d1, d0, t)                        \
do {                                                                        \
    int t##i;                                                               \
    mp_limb_t t##a1 = n1, t##a0 = n0, t##b1 = d1, t##b0 = d0;               \
    mp_limb_t t##r1, t##r0, t##u2, t##u1, t##u0, t##q, t##v1, t##v0;        \
                                                                            \
    FLINT_ASSERT(t##a1 != 0);                                               \
    FLINT_ASSERT(t##b1 != 0);                                               \
    FLINT_ASSERT(t##b1 < t##a1 || (t##b1 == t##a1 && t##b0 < t##a0));       \
                                                                            \
    t##q = 1;                                                               \
                                                                            \
    sub_ddmmss(t##r1,t##r0, t##a1,t##a0, t##b1,t##b0);                      \
                                                                            \
t##subtract:                                                                \
                                                                            \
    for (t##i = 2; t##i <= 4; t##i++)                                       \
    {                                                                       \
        sub_dddmmmsss(t##u2,t##u1,t##u0, 0,t##r1,t##r0, 0,t##b1,t##b0);     \
        if (t##u2 != 0)                                                     \
            goto t##quotient_found;                                         \
        t##q += 1;                                                          \
        t##r0 = t##u0;                                                      \
        t##r1 = t##u1;                                                      \
    }                                                                       \
                                                                            \
    if (t##r1 != 0)                                                         \
    {                                                                       \
        int t##ncnt, t##dcnt;                                               \
        mp_limb_t t##qq = 0;                                                \
                                                                            \
        t##ncnt = flint_clz(t##r1);                                \
        t##dcnt = flint_clz(t##b1);                                \
        t##dcnt -= t##ncnt;                                                 \
        if (t##dcnt <= 0)                                                   \
            goto t##subtract;                                               \
                                                                            \
        t##v1 = (t##b1 << t##dcnt) | (t##b0 >> (FLINT_BITS - t##dcnt));     \
        t##v0 = t##b0 << t##dcnt;                                           \
                                                                            \
        do {                                                                \
            sub_dddmmmsss(t##u2,t##u1,t##u0, 0,t##r1,t##r0, 0,t##v1,t##v0); \
            t##qq = 2*t##qq + 1 + t##u2;                                    \
            t##r1 = t##u2 ? t##r1 : t##u1;                                  \
            t##r0 = t##u2 ? t##r0 : t##u0;                                  \
            t##v0 = (t##v1 << (FLINT_BITS - 1)) | (t##v0 >> 1);             \
            t##v1 = t##v1 >> 1;                                             \
            t##dcnt--;                                                      \
        } while (t##dcnt >= 0);                                             \
                                                                            \
        t##q += t##qq;                                                      \
    }                                                                       \
                                                                            \
t##quotient_found:                                                          \
                                                                            \
    FLINT_ASSERT(t##r1 < t##b1 || (t##r1 == t##b1 && t##r0 < t##b0));       \
                                                                            \
    q0 = t##q;                                                              \
    r0 = t##r0;                                                             \
    r1 = t##r1;                                                             \
} while (0)


#ifdef __cplusplus
}
#endif

#endif
