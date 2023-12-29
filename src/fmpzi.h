/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZI_H
#define FMPZI_H

#ifdef FMPZI_INLINES_C
#define FMPZI_INLINE
#else
#define FMPZI_INLINE static inline
#endif

#include "fmpz.h"

#ifdef __cplusplus
extern "C" {
#endif

#define fmpzi_realref(x) (&((x)->a))
#define fmpzi_imagref(x) (&((x)->b))

FMPZI_INLINE void
fmpzi_init(fmpzi_t x)
{
    fmpz_init(fmpzi_realref(x));
    fmpz_init(fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_clear(fmpzi_t x)
{
    fmpz_clear(fmpzi_realref(x));
    fmpz_clear(fmpzi_imagref(x));
}

FMPZI_INLINE int
fmpzi_equal(const fmpzi_t x, const fmpzi_t y)
{
    return fmpz_equal(fmpzi_realref(x), fmpzi_realref(y)) &&
           fmpz_equal(fmpzi_imagref(x), fmpzi_imagref(y));
}

FMPZI_INLINE void
fmpzi_zero(fmpzi_t x)
{
    fmpz_zero(fmpzi_realref(x));
    fmpz_zero(fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_one(fmpzi_t x)
{
    fmpz_one(fmpzi_realref(x));
    fmpz_zero(fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_set(fmpzi_t res, const fmpzi_t x)
{
    fmpz_set(fmpzi_realref(res), fmpzi_realref(x));
    fmpz_set(fmpzi_imagref(res), fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_conj(fmpzi_t res, const fmpzi_t x)
{
    fmpz_set(fmpzi_realref(res), fmpzi_realref(x));
    fmpz_neg(fmpzi_imagref(res), fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_swap(fmpzi_t x, fmpzi_t y)
{
    FLINT_SWAP(fmpzi_struct, *x, *y);
}

void fmpzi_print(const fmpzi_t x);

FMPZI_INLINE void
fmpzi_set_si_si(fmpzi_t res, slong a, slong b)
{
    fmpz_set_si(fmpzi_realref(res), a);
    fmpz_set_si(fmpzi_imagref(res), b);
}

FMPZI_INLINE void
fmpzi_randtest(fmpzi_t res, flint_rand_t state, mp_bitcnt_t bits)
{
    fmpz_randtest(fmpzi_realref(res), state, bits);
    fmpz_randtest(fmpzi_imagref(res), state, bits);
}

/* Special values */

FMPZI_INLINE int
fmpzi_is_unit(const fmpzi_t x)
{
    if (fmpz_is_zero(fmpzi_imagref(x)))
        return fmpz_is_pm1(fmpzi_realref(x));
    if (fmpz_is_zero(fmpzi_realref(x)))
        return fmpz_is_pm1(fmpzi_imagref(x));
    return 0;
}

FMPZI_INLINE int fmpzi_is_zero(const fmpzi_t x)
{
    return fmpz_is_zero(fmpzi_realref(x)) && fmpz_is_zero(fmpzi_imagref(x));
}

FMPZI_INLINE int fmpzi_is_one(const fmpzi_t x)
{
    return fmpz_is_one(fmpzi_realref(x)) && fmpz_is_zero(fmpzi_imagref(x));
}

/* Norms */

slong fmpzi_bits(const fmpzi_t x);

FMPZI_INLINE void fmpzi_norm(fmpz_t res, const fmpzi_t x)
{
    fmpz_fmma(res, fmpzi_realref(x), fmpzi_realref(x), fmpzi_imagref(x), fmpzi_imagref(x));
}

/* Arithmetic */

FMPZI_INLINE void
fmpzi_neg(fmpzi_t res, const fmpzi_t x)
{
    fmpz_neg(fmpzi_realref(res), fmpzi_realref(x));
    fmpz_neg(fmpzi_imagref(res), fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_add(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
{
    fmpz_add(fmpzi_realref(res), fmpzi_realref(x), fmpzi_realref(y));
    fmpz_add(fmpzi_imagref(res), fmpzi_imagref(x), fmpzi_imagref(y));
}

FMPZI_INLINE void
fmpzi_sub(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
{
    fmpz_sub(fmpzi_realref(res), fmpzi_realref(x), fmpzi_realref(y));
    fmpz_sub(fmpzi_imagref(res), fmpzi_imagref(x), fmpzi_imagref(y));
}

void fmpzi_sqr(fmpzi_t res, const fmpzi_t x);
void fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y);
void fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp);

void fmpzi_mul_i(fmpzi_t z, const fmpzi_t x);
void fmpzi_div_i(fmpzi_t z, const fmpzi_t x);
void fmpzi_mul_i_pow_si(fmpzi_t res, const fmpzi_t z, slong k);

/* todo */
slong fmpzi_canonical_unit_i_pow(const fmpzi_t x);

FMPZI_INLINE void
fmpzi_canonicalise_unit(fmpzi_t res, const fmpzi_t x)
{
    fmpzi_mul_i_pow_si(res, x, fmpzi_canonical_unit_i_pow(x));
}

/* Division */
void fmpzi_divrem(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y);
void fmpzi_divrem_approx(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y);

void fmpzi_divexact(fmpzi_t q, const fmpzi_t x, const fmpzi_t y);

slong fmpzi_remove_one_plus_i(fmpzi_t res, const fmpzi_t x);

/* GCD */
void fmpzi_gcd_euclidean(fmpzi_t res, const fmpzi_t x, const fmpzi_t y);
void fmpzi_gcd_euclidean_improved(fmpzi_t res, const fmpzi_t x, const fmpzi_t y);
void fmpzi_gcd_binary(fmpzi_t res, const fmpzi_t x, const fmpzi_t y);
void fmpzi_gcd_shortest(fmpzi_t g, const fmpzi_t x, const fmpzi_t y);
void fmpzi_gcd(fmpzi_t g, const fmpzi_t x, const fmpzi_t y);

/* Primality check */
int fmpzi_is_prime(const fmpzi_t n);
int fmpzi_is_probabprime(const fmpzi_t n);

#ifdef __cplusplus
}
#endif

#endif
