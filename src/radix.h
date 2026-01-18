/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef RADIX_H
#define RADIX_H

#ifdef RADIX_INLINES_C
#define RADIX_INLINE
#else
#define RADIX_INLINE static inline
#endif

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "nmod.h"
#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    nmod_t b;           /* Digit radix */
    unsigned int exp;   /* B = b^exp */
    nmod_t B;           /* Limb radix */
    /* extended data which may be used by some ring contexts */
    slong trunc_limbs;
    slong trunc_digits;
}
radix_struct;

typedef radix_struct radix_t[1];

#define DIGIT_RADIX(radix) ((radix)->b.n)
#define LIMB_RADIX(radix) ((radix)->B.n)

void radix_init(radix_t radix, ulong b, unsigned int exp);
void radix_clear(radix_t radix);
void radix_init_randtest(radix_t radix, flint_rand_t state);

void radix_rand_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix);
void radix_rand_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix);
void radix_randtest_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix);
void radix_randtest_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix);

ulong radix_neg(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
ulong radix_add(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);
ulong radix_sub(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix);

/* todo
#define radix_add_n(res, a, b, n, radix) radix_add(res, a, n, b, n, radix)
#define radix_sub_n(res, a, b, n, radix) radix_sub(res, a, n, b, n, radix)

RADIX_INLINE ulong
radix_add_1(nn_ptr res, nn_srcptr a, slong n, ulong c, const radix_t radix)
{
    return radix_add(res, a, n, &c, 1, radix);
}

RADIX_INLINE ulong
radix_sub_1(nn_ptr res, nn_srcptr a, slong n, ulong c, const radix_t radix)
{
    return radix_sub(res, a, n, &c, 1, radix);
}
 */

/* Multiplication */

void radix_mulmid_fft_small(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_classical(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_KS(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_naive(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);

RADIX_INLINE void
radix_mulmid(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix)
{
    /* todo: tuning */
    if (bn < 80 || hi - lo < 80)
        radix_mulmid_classical(res, a, an, b, bn, lo, hi, radix);
    else
#if FLINT_HAVE_FFT_SMALL
        radix_mulmid_fft_small(res, a, an, b, bn, lo, hi, radix);
#else
        radix_mulmid_KS(res, a, an, b, bn, lo, hi, radix);
#endif
}

RADIX_INLINE void
radix_mul(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
{
    radix_mulmid(res, a, an, b, bn, 0, an + bn, radix);
}

/* Division */

/* todo: squaring optimisations in all multiplication algorithms */
RADIX_INLINE void
radix_sqr(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    radix_mul(res, a, an, a, an, radix);
}

ulong radix_divrem_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix);
void radix_divexact_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix);

/* Radix conversion */

typedef struct
{
    slong len;
    ulong exps[FLINT_BITS];
    nn_ptr pows[FLINT_BITS];
    slong sizes[FLINT_BITS];
    slong val_limbs[FLINT_BITS];
    nn_ptr buf;
}
radix_powers_struct;

typedef radix_powers_struct radix_powers_t[1];

void radix_powers_clear(radix_powers_t powers);

slong radix_get_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
slong radix_get_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
slong radix_get_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);

slong radix_set_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
slong radix_set_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);
slong radix_set_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix);

slong radix_set_mpn_need_alloc(slong n, const radix_t radix);

/* String conversion */

char * radix_get_str_decimal(char * res, nn_srcptr x, slong n, int negative, const radix_t radix);
char * radix_get_str_sum(char * res, nn_srcptr x, slong n, int negative, int ascending, const radix_t radix);

/* Memory-managed integers */

typedef struct
{
    nn_ptr d;
    slong alloc;
    slong size;
}
radix_integer_struct;

typedef radix_integer_struct radix_integer_t[1];

void gr_ctx_init_radix_integer(gr_ctx_t ctx, ulong b, unsigned int exp);

void radix_integer_init(radix_integer_t res, const radix_t radix);
void radix_integer_clear(radix_integer_t res, const radix_t radix);
nn_ptr radix_integer_fit_limbs(radix_integer_t res, slong nlimbs, const radix_t radix);
void radix_integer_zero(radix_integer_t res, const radix_t radix);
void radix_integer_randtest_limbs(radix_integer_t res, flint_rand_t state, slong max_limbs, const radix_t radix);
void radix_integer_one(radix_integer_t res, const radix_t radix);
void radix_integer_neg_one(radix_integer_t res, const radix_t radix);
int radix_integer_is_zero(const radix_integer_t x, const radix_t radix);
int radix_integer_is_one(const radix_integer_t x, const radix_t radix);
int radix_integer_is_neg_one(const radix_integer_t x, const radix_t radix);
int radix_integer_equal(const radix_integer_t x, const radix_integer_t y, const radix_t radix);
int radix_integer_cmp(const radix_integer_t x, const radix_integer_t y, const radix_t radix);
int radix_integer_cmpabs(const radix_integer_t x, const radix_integer_t y, const radix_t radix);
void radix_integer_set(radix_integer_t res, const radix_integer_t x, const radix_t radix);
void radix_integer_set_ui(radix_integer_t res, ulong x, const radix_t radix);
void radix_integer_set_si(radix_integer_t res, slong x, const radix_t radix);
void radix_integer_set_fmpz(radix_integer_t res, const fmpz_t x, const radix_t radix);
void radix_integer_get_fmpz(fmpz_t res, const radix_integer_t x, const radix_t radix);
void radix_integer_neg(radix_integer_t res, const radix_integer_t x, const radix_t radix);
void radix_integer_abs(radix_integer_t res, const radix_integer_t x, const radix_t radix);
int radix_integer_sgn(const radix_integer_t x, const radix_t radix);
void radix_integer_add(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix);
void radix_integer_sub(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix);
void radix_integer_mul(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix);

#ifdef __cplusplus
}
#endif

#endif
