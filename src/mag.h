/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MAG_H
#define MAG_H

#ifdef MAG_INLINES_C
#define MAG_INLINE
#else
#define MAG_INLINE static inline
#endif

#include "fmpz_extras.h"
#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define LIMB_ONE ((mp_limb_t) 1)
#define LIMB_ONES (-(mp_limb_t) 1)
#define LIMB_TOP (((mp_limb_t) 1) << (FLINT_BITS - 1))
#define MASK_LIMB(n, c) ((n) & (LIMB_ONES << (c)))

#define MAG_MAX_LAGOM_EXP (COEFF_MAX / 4)
#define MAG_MIN_LAGOM_EXP (-MAG_MAX_LAGOM_EXP)

#define ADD2_FAST_MAX (COEFF_MAX / 4)
#define ADD2_FAST_MIN (-ADD2_FAST_MAX)

/* TODO: rename these and move them to fmpz_extras */

static inline void
_fmpz_set_fast(fmpz_t f, const fmpz_t g)
{
    if (!COEFF_IS_MPZ(*f) && !COEFF_IS_MPZ(*g))
        *f = *g;
    else
        fmpz_set(f, g);
}

static inline void
_fmpz_add_fast(fmpz_t z, const fmpz_t x, slong c)
{
    fmpz ze, xe;

    ze = *z;
    xe = *x;

    if (!COEFF_IS_MPZ(ze) && (xe > ADD2_FAST_MIN && xe < ADD2_FAST_MAX))
        *z = xe + c;
    else
        fmpz_add_si(z, x, c);
}

static inline void
_fmpz_add2_fast(fmpz_t z, const fmpz_t x, const fmpz_t y, slong c)
{
    fmpz ze, xe, ye;

    ze = *z;
    xe = *x;
    ye = *y;

    if (!COEFF_IS_MPZ(ze) && (xe > ADD2_FAST_MIN && xe < ADD2_FAST_MAX) &&
                             (ye > ADD2_FAST_MIN && ye < ADD2_FAST_MAX))
    {
        *z = xe + ye + c;
    }
    else
    {
        fmpz_add(z, x, y);
        fmpz_add_si(z, z, c);
    }
}

static inline void
_fmpz_sub2_fast(fmpz_t z, const fmpz_t x, const fmpz_t y, slong c)
{
    fmpz ze, xe, ye;

    ze = *z;
    xe = *x;
    ye = *y;

    if (!COEFF_IS_MPZ(ze) && (xe > ADD2_FAST_MIN && xe < ADD2_FAST_MAX) &&
                             (ye > ADD2_FAST_MIN && ye < ADD2_FAST_MAX))
    {
        *z = xe - ye + c;
    }
    else
    {
        fmpz_sub(z, x, y);
        fmpz_add_si(z, z, c);
    }
}


#define MAG_EXP_POS_INF (COEFF_MIN+1)

/* Finite and with lagom big exponents. */
#define MAG_IS_LAGOM(x) (MAG_EXP(x) >= MAG_MIN_LAGOM_EXP && \
                         MAG_EXP(x) <= MAG_MAX_LAGOM_EXP)

#define MAG_EXPREF(x) (&(x)->exp)
#define MAG_EXP(x) ((x)->exp)
#define MAG_MAN(x) ((x)->man)

#define MAG_BITS 30

#define MAG_ONE_HALF (UWORD(1) << (MAG_BITS - 1))

static inline mp_limb_t
__mag_fixmul32(mp_limb_t x, mp_limb_t y)
{
    mp_limb_t u, v;
    umul_ppmm(u, v, x, y);
    return (u << (32 - MAG_BITS)) | (v >> MAG_BITS);
}

#if FLINT_BITS == 64
#define MAG_FIXMUL(x, y) (((x) * (y)) >> MAG_BITS)
#else
#define MAG_FIXMUL(x, y) __mag_fixmul32((x), (y))
#endif

#define MAG_CHECK_BITS(rr) \
    if (MAG_MAN(rr) != 0 && FLINT_BIT_COUNT(MAG_MAN(rr)) != MAG_BITS) \
    { \
        flint_throw(FLINT_ERROR, "FAIL: wrong number of bits in mantissa!\n"); \
    }

/* Note: assumes mantissa either has the right number of bits, or
   one more bit (but in that case must not be all ones, as that
   would round up to one extra bit again, requiring a second
   correction). */

#define MAG_ADJUST_ONE_TOO_LARGE(x) \
    do { \
        mp_limb_t __t = MAG_MAN(x) >> MAG_BITS; \
        MAG_MAN(x) = (MAG_MAN(x) >> __t) + (__t & MAG_MAN(x)); \
        if (__t) \
            fmpz_add_ui(MAG_EXPREF(x), MAG_EXPREF(x), __t); \
    } while (0)

#define MAG_FAST_ADJUST_ONE_TOO_LARGE(x) \
    do { \
        mp_limb_t __t = MAG_MAN(x) >> MAG_BITS; \
        MAG_MAN(x) = (MAG_MAN(x) >> __t) + (__t & MAG_MAN(x)); \
        MAG_EXP(x) += __t; \
    } while (0)

#define MAG_ADJUST_ONE_TOO_SMALL(x) \
    do { \
        mp_limb_t __t = !(MAG_MAN(x) >> (MAG_BITS - 1)); \
        MAG_MAN(x) = (MAG_MAN(x) << __t); \
        if (__t) \
            fmpz_sub_ui(MAG_EXPREF(x), MAG_EXPREF(x), __t); \
    } while (0)

#define MAG_FAST_ADJUST_ONE_TOO_SMALL(x) \
    do { \
        mp_limb_t __t = !(MAG_MAN(x) >> (MAG_BITS - 1)); \
        MAG_MAN(x) = (MAG_MAN(x) << __t); \
        MAG_EXP(x) -= __t; \
    } while (0)


MAG_INLINE void
mag_init(mag_t x)
{
    fmpz_init(MAG_EXPREF(x));
    MAG_MAN(x) = 0;
}

MAG_INLINE void
mag_init_set(mag_t x, const mag_t y)
{
    fmpz_init_set(MAG_EXPREF(x), MAG_EXPREF(y));
    MAG_MAN(x) = MAG_MAN(y);
}

void mag_clear(mag_t x);

MAG_INLINE void
mag_swap(mag_t x, mag_t y)
{
    mag_struct t = *x;
    *x = *y;
    *y = t;
}

MAG_INLINE void
mag_set(mag_t x, const mag_t y)
{
    _fmpz_set_fast(MAG_EXPREF(x), MAG_EXPREF(y));
    x->man = y->man;
}

MAG_INLINE void
mag_zero(mag_t x)
{
    fmpz_zero(MAG_EXPREF(x));
    MAG_MAN(x) = 0;
}

MAG_INLINE void
mag_one(mag_t x)
{
    fmpz_one(MAG_EXPREF(x));
    MAG_MAN(x) = MAG_ONE_HALF;
}

MAG_INLINE int
mag_is_special(const mag_t x)
{
    return MAG_MAN(x) == 0;
}

MAG_INLINE int
mag_is_zero(const mag_t x)
{
    return (MAG_MAN(x) == 0) && (MAG_EXP(x) == 0);
}

MAG_INLINE void
mag_inf(mag_t x)
{
    fmpz_clear(MAG_EXPREF(x));
    MAG_EXP(x) = MAG_EXP_POS_INF;
    MAG_MAN(x) = 0;
}

MAG_INLINE int
mag_is_inf(const mag_t x)
{
    return (MAG_MAN(x) == 0) && (MAG_EXP(x) != 0);
}

MAG_INLINE int
mag_is_finite(const mag_t x)
{
    return !mag_is_inf(x);
}

MAG_INLINE int
mag_equal(const mag_t x, const mag_t y)
{
    return (MAG_MAN(x) == MAG_MAN(y))
        && fmpz_equal(MAG_EXPREF(x), MAG_EXPREF(y));
}

/* general versions */

void mag_mul(mag_t z, const mag_t x, const mag_t y);

void mag_mul_lower(mag_t z, const mag_t x, const mag_t y);

void mag_addmul(mag_t z, const mag_t x, const mag_t y);

void mag_add_2exp_fmpz(mag_t z, const mag_t x, const fmpz_t e);

void mag_add(mag_t z, const mag_t x, const mag_t y);

void mag_add_lower(mag_t z, const mag_t x, const mag_t y);

void mag_add_ui(mag_t z, const mag_t x, ulong y);
void mag_add_ui_lower(mag_t res, const mag_t x, ulong y);

void mag_add_ui_2exp_si(mag_t z, const mag_t x, ulong y, slong e);

void mag_div(mag_t z, const mag_t x, const mag_t y);
void mag_div_lower(mag_t z, const mag_t x, const mag_t y);

MAG_INLINE void
mag_inv(mag_t res, const mag_t x)
{
    mag_t t;
    mag_init(t); /* no need to free */
    mag_one(t);
    mag_div(res, t, x);
}

MAG_INLINE void
mag_inv_lower(mag_t res, const mag_t x)
{
    mag_t t;
    mag_init(t); /* no need to free */
    mag_one(t);
    mag_div_lower(res, t, x);
}

void mag_mul_2exp_si(mag_t z, const mag_t x, slong y);

void mag_mul_2exp_fmpz(mag_t z, const mag_t x, const fmpz_t y);

void mag_sub(mag_t z, const mag_t x, const mag_t y);

void mag_sub_lower(mag_t z, const mag_t x, const mag_t y);

/* Fast versions (no infs/nans, small exponents). Note that this
   applies to outputs too! */

MAG_INLINE void
mag_fast_init_set(mag_t x, const mag_t y)
{
    MAG_EXP(x) = MAG_EXP(y);
    MAG_MAN(x) = MAG_MAN(y);
}

MAG_INLINE void
mag_fast_zero(mag_t x)
{
    MAG_EXP(x) = 0;
    MAG_MAN(x) = 0;
}

MAG_INLINE int
mag_fast_is_zero(const mag_t x)
{
    return MAG_MAN(x) == 0;
}

MAG_INLINE void
mag_fast_mul(mag_t z, const mag_t x, const mag_t y)
{
    if (MAG_MAN(x) == 0 || MAG_MAN(y) == 0)
    {
        mag_fast_zero(z);
    }
    else
    {
        MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + LIMB_ONE;
        MAG_EXP(z) = MAG_EXP(x) + MAG_EXP(y);
        MAG_FAST_ADJUST_ONE_TOO_SMALL(z);
    }
}

MAG_INLINE void
mag_fast_mul_2exp_si(mag_t z, const mag_t x, slong y)
{
    if (MAG_MAN(x) == 0)
    {
        mag_fast_zero(z);
    }
    else
    {
        MAG_MAN(z) = MAG_MAN(x);
        MAG_EXP(z) = MAG_EXP(x) + y;
    }
}

MAG_INLINE void
mag_fast_addmul(mag_t z, const mag_t x, const mag_t y)
{
    if (MAG_MAN(z) == 0)
    {
        mag_fast_mul(z, x, y);
    }
    else if (MAG_MAN(x) == 0 || MAG_MAN(y) == 0)
    {
        return;
    }
    else
    {
        slong shift, e;

        /* x*y < 2^e */
        e = MAG_EXP(x) + MAG_EXP(y);
        shift = MAG_EXP(z) - e;

        if (shift >= 0)
        {
            if (shift >= MAG_BITS)
                MAG_MAN(z)++;
            else
                MAG_MAN(z) = MAG_MAN(z) + (MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) >> shift) + 1;
        }
        else
        {
            shift = -shift;
            MAG_EXP(z) = e;

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + 2;
            else
                MAG_MAN(z) = MAG_FIXMUL(MAG_MAN(x), MAG_MAN(y)) + (MAG_MAN(z) >> shift) + 2;

            MAG_FAST_ADJUST_ONE_TOO_SMALL(z);
        }

        MAG_FAST_ADJUST_ONE_TOO_LARGE(z);
    }
}

MAG_INLINE void
mag_fast_add_2exp_si(mag_t z, const mag_t x, slong e)
{
    /* Must be zero */
    if (mag_is_special(x))
    {
        MAG_MAN(z) = MAG_ONE_HALF;
        MAG_EXP(z) = e + 1;
    }
    else
    {
        slong shift;
        shift = MAG_EXP(x) - e;

        if (shift > 0)
        {
            MAG_EXP(z) = MAG_EXP(x);

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_MAN(x) + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_MAN(x) + (LIMB_ONE << (MAG_BITS - shift));
        }
        else
        {
            shift = -shift;

            MAG_EXP(z) = e + 1;

            if (shift >= MAG_BITS)
                MAG_MAN(z) = MAG_ONE_HALF + LIMB_ONE;
            else
                MAG_MAN(z) = MAG_ONE_HALF + (MAG_MAN(x) >> (shift + 1)) + LIMB_ONE;
        }

        MAG_FAST_ADJUST_ONE_TOO_LARGE(z);
    }
}

/* requires that x is positive and finite */
#define MAG_SET_D_2EXP(man, exp, x, xexp) \
    do { \
        int __cexp; \
        double __x; \
        int __fix; \
        mp_limb_t __man; \
        __x = frexp((x), &__cexp); \
        __man = (mp_limb_t)(__x * (double)(LIMB_ONE << MAG_BITS)) + 1; \
        __fix = __man >> (MAG_BITS); \
        __man = (__man >> __fix) + __fix; \
        (man) = __man; \
        (exp) = (xexp) + __cexp + __fix; \
    } while (0);

/* requires that x is positive and finite */
#define MAG_SET_D_2EXP_LOWER(man, exp, x, xexp) \
    do { \
        int __cexp; \
        double __x; \
        int __fix; \
        mp_limb_t __man; \
        __x = frexp((x), &__cexp); \
        __man = (mp_limb_t)(__x * (double)(LIMB_ONE << MAG_BITS)) - 1; \
        __fix = __man < MAG_ONE_HALF; \
        __man = (__man << __fix); \
        (man) = __man; \
        (exp) = (xexp) + __cexp - __fix; \
    } while (0);

void mag_set_d(mag_t z, double x);
void mag_set_d_lower(mag_t z, double x);
void mag_set_d_2exp_fmpz(mag_t z, double c, const fmpz_t exp);
void mag_set_d_2exp_fmpz_lower(mag_t z, double c, const fmpz_t exp);

void mag_set_fmpz_2exp_fmpz(mag_t z, const fmpz_t man, const fmpz_t exp);

void mag_randtest_special(mag_t x, flint_rand_t state, slong expbits);

void mag_randtest(mag_t x, flint_rand_t state, slong expbits);

#ifdef FLINT_HAVE_FILE
void mag_fprint(FILE * file, const mag_t x);
void mag_fprintd(FILE * file, const mag_t x, slong d);
#endif

void mag_print(const mag_t x);
void mag_printd(const mag_t x, slong d);

void mag_get_fmpq(fmpq_t y, const mag_t x);

void mag_get_fmpz(fmpz_t res, const mag_t x);
void mag_get_fmpz_lower(fmpz_t res, const mag_t x);

int mag_cmp(const mag_t x, const mag_t y);

int mag_cmp_2exp_si(const mag_t x, slong e);

MAG_INLINE void
mag_min(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_cmp(x, y) <= 0)
        mag_set(z, x);
    else
        mag_set(z, y);
}

MAG_INLINE void
mag_max(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_cmp(x, y) >= 0)
        mag_set(z, x);
    else
        mag_set(z, y);
}

MAG_INLINE mag_ptr
_mag_vec_init(slong n)
{
    slong i;
    mag_ptr v = (mag_ptr) flint_malloc(sizeof(mag_struct) * n);

    for (i = 0; i < n; i++)
        mag_init(v + i);

    return v;
}

MAG_INLINE void
_mag_vec_clear(mag_ptr v, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        mag_clear(v + i);
    flint_free(v);
}

double mag_get_d(const mag_t z);

double mag_get_d_log2_approx(const mag_t x);

/* TODO: document */
double mag_d_log_upper_bound(double x);
double mag_d_log_lower_bound(double x);

void mag_log1p(mag_t z, const mag_t x);

void mag_log_ui(mag_t t, ulong n);

void mag_log(mag_t z, const mag_t x);
void mag_log_lower(mag_t z, const mag_t x);
void mag_neg_log(mag_t z, const mag_t x);
void mag_neg_log_lower(mag_t z, const mag_t x);

void mag_exp(mag_t y, const mag_t x);
void mag_exp_lower(mag_t y, const mag_t x);

void mag_expinv(mag_t res, const mag_t x);
void mag_expinv_lower(mag_t y, const mag_t x);

void mag_expm1(mag_t y, const mag_t x);
void mag_exp_tail(mag_t z, const mag_t x, ulong N);

void mag_sinh(mag_t y, const mag_t x);
void mag_sinh_lower(mag_t y, const mag_t x);

void mag_cosh(mag_t y, const mag_t x);
void mag_cosh_lower(mag_t y, const mag_t x);

void mag_pow_ui(mag_t z, const mag_t x, ulong e);
void mag_pow_ui_lower(mag_t z, const mag_t x, ulong e);
void mag_pow_fmpz(mag_t z, const mag_t x, const fmpz_t e);
void mag_pow_fmpz_lower(mag_t z, const mag_t x, const fmpz_t e);

void mag_const_pi(mag_t res);
void mag_const_pi_lower(mag_t res);

void mag_atan(mag_t res, const mag_t x);
void mag_atan_lower(mag_t res, const mag_t x);

void mag_fac_ui(mag_t z, ulong n);
void mag_rfac_ui(mag_t z, ulong n);
void mag_bin_uiui(mag_t res, ulong n, ulong k);

/* TODO: test */
void mag_bernoulli_div_fac_ui(mag_t z, ulong n);

/* TODO: test */
void mag_set_fmpz_2exp_fmpz_lower(mag_t z, const fmpz_t man, const fmpz_t exp);

void mag_sqrt(mag_t y, const mag_t x);
void mag_sqrt_lower(mag_t y, const mag_t x);
void mag_rsqrt(mag_t y, const mag_t x);
void mag_rsqrt_lower(mag_t y, const mag_t x);

void mag_root(mag_t y, const mag_t x, ulong n);

void mag_hypot(mag_t z, const mag_t x, const mag_t y);

void mag_binpow_uiui(mag_t b, ulong m, ulong n);

void mag_polylog_tail(mag_t u, const mag_t z, slong sigma, ulong d, ulong N);

void mag_geom_series(mag_t res, const mag_t x, ulong n);

void mag_hurwitz_zeta_uiui(mag_t res, ulong s, ulong a);

void mag_set_ui(mag_t z, ulong x);
void mag_set_ui_lower(mag_t z, ulong x);

/* TODO: test functions below */
void mag_set_ui_2exp_si(mag_t z, ulong v, slong e);

MAG_INLINE void
mag_set_fmpz(mag_t z, const fmpz_t x)
{
    fmpz_t exp;
    *exp = 0;
    mag_set_fmpz_2exp_fmpz(z, x, exp);
}

MAG_INLINE void
mag_set_fmpz_lower(mag_t z, const fmpz_t x)
{
    fmpz_t exp;
    *exp = 0;
    mag_set_fmpz_2exp_fmpz_lower(z, x, exp);
}

MAG_INLINE void
mag_mul_ui(mag_t z, const mag_t x, ulong y)
{
    mag_t t;
    mag_init(t);  /* no need to clear */
    mag_set_ui(t, y);
    mag_mul(z, x, t);
}

MAG_INLINE void
mag_mul_ui_lower(mag_t z, const mag_t x, ulong y)
{
    mag_t t;
    mag_init(t);  /* no need to clear */
    mag_set_ui_lower(t, y);
    mag_mul_lower(z, x, t);
}

MAG_INLINE void
mag_mul_fmpz(mag_t z, const mag_t x, const fmpz_t y)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz(t, y);
    mag_mul(z, x, t);
    mag_clear(t);
}

MAG_INLINE void
mag_mul_fmpz_lower(mag_t z, const mag_t x, const fmpz_t y)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz_lower(t, y);
    mag_mul_lower(z, x, t);
    mag_clear(t);
}

MAG_INLINE void
mag_div_ui(mag_t z, const mag_t x, ulong y)
{
    mag_t t;
    mag_init(t);  /* no need to clear */
    mag_set_ui_lower(t, y);
    mag_div(z, x, t);
}

MAG_INLINE void
mag_div_fmpz(mag_t z, const mag_t x, const fmpz_t y)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz_lower(t, y);
    mag_div(z, x, t);
    mag_clear(t);
}

MAG_INLINE slong mag_allocated_bytes(const mag_t x)
{
    return fmpz_allocated_bytes(MAG_EXPREF(x));
}

int mag_load_str(mag_t res, const char * data);
char * mag_dump_str(const mag_t x);

#ifdef FLINT_HAVE_FILE
int mag_load_file(mag_t res, FILE * stream);
int mag_dump_file(FILE * stream, const mag_t x);
#endif

#ifdef __cplusplus
}
#endif

#endif

