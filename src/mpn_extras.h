/*
    Copyright (C) 2010 Fredrik Johansson

    2x2 mul code taken from MPFR 2.3.0
    (Copyright (C) 1991-2007 Free Software Foundation, Inc.)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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

/* Macros for few-limb arithmetic */

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

/* {s0,s1,s2} = u[0]v[n-1] + u[1]v[n-2] + ... */
/* Assumes n >= 2 */
#define NN_DOTREV_S3_1X1(s2, s1, s0, u, v, n) \
    do { \
        mp_limb_t __dt0, __dt1, __ds0, __ds1, __ds2; \
        slong __i; \
        FLINT_ASSERT((n) >= 2); \
        umul_ppmm(__ds1, __ds0, (u)[0], (v)[(n) - 1]); \
        umul_ppmm(__dt1, __dt0, (u)[1], (v)[(n) - 2]); \
        add_sssaaaaaa(__ds2, __ds1, __ds0, 0, __ds1, __ds0, 0, __dt1, __dt0); \
        for (__i = 2; __i < (n); __i++) \
        { \
            umul_ppmm(__dt1, __dt0, (u)[i], (v)[(n) - 1 - __i]); \
            add_sssaaaaaa(__ds2, __ds1, __ds0, __ds2, __ds1, __ds0, 0, __dt1, __dt0); \
        } \
        (s0) = __ds0; (s1) = __ds1; (s2) = __ds2; \
    } while (0) \

/* {r0,r1,r2} = {s0,s1,s2} + u[0]v[n-1] + u[1]v[n-2] + ... */
/* Assumes n >= 1. May have s2 != 0, but the final sum is assumed to fit in 3 limbs. */
#define NN_DOTREV_S3_A3_1X1(r2, r1, r0, s2, s1, s0, u, v, n) \
    do { \
        mp_limb_t __dt0, __dt1, __ds0, __ds1, __ds2; \
        slong __i; \
        FLINT_ASSERT((n) >= 1); \
        __ds0 = (s0); __ds1 = (s1); __ds2 = (s2); \
        for (__i = 0; __i < (n); __i++) \
        { \
            umul_ppmm(__dt1, __dt0, (u)[__i], (v)[(n) - 1 - __i]); \
            add_sssaaaaaa(__ds2, __ds1, __ds0, __ds2, __ds1, __ds0, 0, __dt1, __dt0); \
        } \
        (r0) = __ds0; (r1) = __ds1; (r2) = __ds2; \
    } while (0) \

#define NN_MUL_1X1 umul_ppmm

/* {r0,r1} = {s0,s1} + x * y, with no carry-out. */
#define NN_ADDMUL_S2_A2_1X1(r1, r0, s1, s0, x, y) \
    do { \
        mp_limb_t __dt0, __dt1; \
        umul_ppmm(__dt1, __dt0, (x), (y)); \
        add_ssaaaa(r1, r0, s1, s0, __dt1, __dt0); \
    } while (0); \

double
flint_mpn_get_d(mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp);

/* General multiplication ****************************************************/

#ifdef FLINT_HAVE_FFT_SMALL
# define FLINT_FFT_MUL_THRESHOLD FLINT_FFT_SMALL_THRESHOLD
# define FLINT_FFT_SQR_THRESHOLD (2 * FLINT_FFT_SMALL_THRESHOLD)
#else
/* FLINT's FFT can beat GMP below this threshold but apparently
   not consistently. Something needs retuning? */
# define FLINT_FFT_MUL_THRESHOLD 32000
# define FLINT_FFT_SQR_THRESHOLD 32000
#endif

#if FLINT_HAVE_ADX
# define FLINT_MPN_MUL_FUNC_TAB_WIDTH 17
# define FLINT_MPN_SQR_FUNC_TAB_WIDTH 14
# define FLINT_HAVE_MUL_FUNC(n, m) ((n) <= 16)
# define FLINT_HAVE_MUL_N_FUNC(n) ((n) <= 16)
# define FLINT_HAVE_SQR_FUNC(n) ((n) <= FLINT_MPN_SQR_FUNC_TAB_WIDTH)
#else
# define FLINT_MPN_MUL_FUNC_TAB_WIDTH 8
# define FLINT_HAVE_MUL_FUNC(n, m) ((n) <= 7 || ((n) <= 14 && (m) == 1))
# define FLINT_HAVE_MUL_N_FUNC(n) ((n) <= 7)
# define FLINT_HAVE_SQR_FUNC(n) (0)
#endif

#define FLINT_MUL_USE_FUNC_TAB 1

typedef mp_limb_t (* flint_mpn_mul_func_t)(mp_ptr, mp_srcptr, mp_srcptr);
typedef mp_limb_t (* flint_mpn_sqr_func_t)(mp_ptr, mp_srcptr);

FLINT_DLL extern const flint_mpn_mul_func_t flint_mpn_mul_func_tab[][FLINT_MPN_MUL_FUNC_TAB_WIDTH];
FLINT_DLL extern const flint_mpn_mul_func_t flint_mpn_mul_n_func_tab[];
FLINT_DLL extern const flint_mpn_sqr_func_t flint_mpn_sqr_func_tab[];

mp_limb_t flint_mpn_mul_basecase(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t xn, mp_size_t yn);
mp_limb_t flint_mpn_sqr_basecase(mp_ptr r, mp_srcptr x, mp_size_t n);

void flint_mpn_mul_toom22(mp_ptr pp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn, mp_ptr scratch);

mp_limb_t _flint_mpn_mul(mp_ptr r, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn);
void _flint_mpn_mul_n(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n);
mp_limb_t _flint_mpn_sqr(mp_ptr r, mp_srcptr x, mp_size_t n);

MPN_EXTRAS_INLINE mp_limb_t
flint_mpn_mul(mp_ptr r, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)
{
    FLINT_ASSERT(xn >= yn);
    FLINT_ASSERT(yn >= 1);
    FLINT_ASSERT(r != x);
    FLINT_ASSERT(r != y);

    if (FLINT_MUL_USE_FUNC_TAB && FLINT_HAVE_MUL_FUNC(xn, yn))
        return flint_mpn_mul_func_tab[xn][yn](r, x, y);
    else
        return _flint_mpn_mul(r, x, xn, y, yn);
}

MPN_EXTRAS_INLINE void
flint_mpn_mul_n(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r != x);
    FLINT_ASSERT(r != y);

    if (FLINT_MUL_USE_FUNC_TAB && FLINT_HAVE_MUL_N_FUNC(n))
        flint_mpn_mul_n_func_tab[n](r, x, y);
    else
        _flint_mpn_mul_n(r, x, y, n);
}

MPN_EXTRAS_INLINE mp_limb_t
flint_mpn_sqr(mp_ptr r, mp_srcptr x, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_MUL_USE_FUNC_TAB && FLINT_HAVE_SQR_FUNC(n))
    {
        /* NOTE: Aliasing allowed */
        return flint_mpn_sqr_func_tab[n](r, x);
    }
    else
    {
        FLINT_ASSERT(r != x);
        return _flint_mpn_sqr(r, x, n);
    }
}

/* Like flint_mpn_mul but allow operands in either order, completely
   inline some small products, and also check for squaring. */
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
        flint_mpn_mul((_z), (_x), (_xn), (_y), (_yn)); \
    } \
    else \
    { \
        flint_mpn_mul((_z), (_y), (_yn), (_x), (_xn)); \
    }

/* High multiplication *******************************************************/

#define FLINT_HAVE_MULHIGH_FUNC(n) ((n) <= FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH)
#define FLINT_HAVE_SQRHIGH_FUNC(n) ((n) <= FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH)
#define FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n) ((n) <= FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH)

struct mp_limb_pair_t { mp_limb_t m1; mp_limb_t m2; };
typedef struct mp_limb_pair_t (* flint_mpn_mulhigh_normalised_func_t)(mp_ptr, mp_srcptr, mp_srcptr);

FLINT_DLL extern const flint_mpn_mul_func_t flint_mpn_mulhigh_func_tab[];
FLINT_DLL extern const flint_mpn_sqr_func_t flint_mpn_sqrhigh_func_tab[];
FLINT_DLL extern const flint_mpn_mulhigh_normalised_func_t flint_mpn_mulhigh_normalised_func_tab[];

#if FLINT_HAVE_ADX
# define FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH 12
# define FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH 8
# define FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH 12

# define FLINT_HAVE_NATIVE_MPN_MULHIGH_BASECASE 1
# define FLINT_HAVE_NATIVE_MPN_SQRHIGH_BASECASE 1

/* NOTE: This function only works for n >= 6 */
mp_limb_t _flint_mpn_mulhigh_basecase(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

/* NOTE: These two functions only works for n >= 8 */
mp_limb_t _flint_mpn_sqrhigh_basecase_even(mp_ptr, mp_srcptr, mp_size_t);
mp_limb_t _flint_mpn_sqrhigh_basecase_odd(mp_ptr, mp_srcptr, mp_size_t);

/* TODO: Proceed with higher cases */
MPN_EXTRAS_INLINE
mp_limb_t flint_mpn_mulhigh_basecase(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_FUNC(n)) /* NOTE: Aliasing allowed here */
        return flint_mpn_mulhigh_func_tab[n](rp, xp, yp);
    else
    {
        FLINT_ASSERT(rp != xp && rp != yp);
        return _flint_mpn_mulhigh_basecase(rp, xp, yp, n);
    }
}

/* TODO: Proceed with higher cases */
MPN_EXTRAS_INLINE
mp_limb_t flint_mpn_sqrhigh_basecase(mp_ptr rp, mp_srcptr xp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_SQRHIGH_FUNC(n)) /* NOTE: Aliasing allowed here */
        return flint_mpn_sqrhigh_func_tab[n](rp, xp);
    else
    {
        FLINT_ASSERT(rp != xp);
        if (n & 1)
            return _flint_mpn_sqrhigh_basecase_odd(rp, xp, n >> 1);
        else
            return _flint_mpn_sqrhigh_basecase_even(rp, xp, n >> 1);
    }
}

MPN_EXTRAS_INLINE
struct mp_limb_pair_t flint_mpn_mulhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n >= 1);

    if (FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n))
        return flint_mpn_mulhigh_normalised_func_tab[n](rp, xp, yp);
    else
    {
        struct mp_limb_pair_t ret;

        FLINT_ASSERT(rp != xp && rp != yp);

        /* TODO */
        /* ret.m1 = flint_mpn_mulhigh(rp, xp, yp, n); */
        ret.m1 = flint_mpn_mulhigh_basecase(rp, xp, yp, n);

        if (rp[n - 1] >> (FLINT_BITS - 1))
        {
            ret.m2 = 0;
        }
        else
        {
            ret.m2 = 1;
            mpn_lshift(rp, rp, n, 1);
            rp[0] |= (ret.m1 >> (FLINT_BITS - 1));
            ret.m1 <<= 1;
        }

        return ret;
    }
}
#else
# define FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH 0
# define FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH 0
# define FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH 0
#endif

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

mp_limb_t flint_mpn_sumdiff_n(mp_ptr s, mp_ptr d, mp_srcptr x, mp_srcptr y, mp_size_t n);

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
