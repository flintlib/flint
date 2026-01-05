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

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    nmod_t b;           /* Digit radix */
    unsigned int exp;   /* B = b^exp */
    nmod_t B;           /* Limb radix */
}
radix_struct;

typedef radix_struct radix_t[1];

#define DIGIT_RADIX(radix) ((radix)->b.n)
#define LIMB_RADIX(radix) ((radix)->B.n)

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

void radix_mulmid_fft_small(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_classical(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_KS(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);
void radix_mulmid_naive(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix);

RADIX_INLINE void
radix_mulmid(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix)
{
    /* todo: tuning */
    if (bn < 30 || hi - lo < 30)
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

/* todo: squaring optimisations in all multiplication algorithms */
RADIX_INLINE void
radix_sqr(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
{
    radix_mul(res, a, an, a, an, radix);
}

ulong radix_divrem_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix);
void radix_divexact_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix);

#ifdef __cplusplus
}
#endif

#endif
