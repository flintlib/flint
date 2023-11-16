/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_EXTRAS_H
#define FMPZ_EXTRAS_H

#include "fmpz.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline void
fmpz_add_inline(fmpz_t z, const fmpz_t x, const fmpz_t y)
{
    fmpz f, g;

    f = *x;
    g = *y;

    if (!COEFF_IS_MPZ(f) && !COEFF_IS_MPZ(g))
        fmpz_set_si(z, f + g);
    else
        fmpz_add(z, x, y);
}

static inline void
fmpz_add_si_inline(fmpz_t z, const fmpz_t x, slong y)
{
    fmpz f;

    f = *x;

    if (!COEFF_IS_MPZ(f) && (COEFF_MIN <= y && y <= COEFF_MAX))
        fmpz_set_si(z, f + y);
    else
        fmpz_add_si(z, x, y);
}

static inline void
fmpz_sub_si_inline(fmpz_t z, const fmpz_t x, slong y)
{
    fmpz f;

    f = *x;

    if (!COEFF_IS_MPZ(f) && (COEFF_MIN <= y && y <= COEFF_MAX))
        fmpz_set_si(z, f - y);
    else
        fmpz_sub_si(z, x, y);
}

static inline void
fmpz_add_ui_inline(fmpz_t z, const fmpz_t x, ulong y)
{
    fmpz f = *x;

    if (!COEFF_IS_MPZ(f) && y <= COEFF_MAX)
        fmpz_set_si(z, f + y);
    else
        fmpz_add_ui(z, x, y);
}

static inline void
fmpz_add2_fmpz_si_inline(fmpz_t z, const fmpz_t x, const fmpz_t y, slong c)
{
    fmpz f, g, h;

    f = *x;
    g = *y;

    if (!COEFF_IS_MPZ(f) && !COEFF_IS_MPZ(g))
    {
        h = f + g;

        if ((COEFF_MIN <= h && h <= COEFF_MAX) &&
            (COEFF_MIN <= c && c <= COEFF_MAX))
        {
            fmpz_set_si(z, h + c);
            return;
        }
    }

    fmpz_add(z, x, y);
    fmpz_add_si(z, z, c);
}

static inline void
fmpz_set_mpn_large(fmpz_t z, mp_srcptr src, mp_size_t n, int negative)
{
    __mpz_struct * zz;
    zz = _fmpz_promote(z);

    if (zz->_mp_alloc < n)
        mpz_realloc2(zz, n * FLINT_BITS);

    flint_mpn_copyi(zz->_mp_d, src, n);
    zz->_mp_size = negative ? -n : n;
}

static inline void
fmpz_adiv_q_2exp(fmpz_t z, const fmpz_t x, flint_bitcnt_t exp)
{
    int sign = fmpz_sgn(x);

    if (sign > 0)
        fmpz_cdiv_q_2exp(z, x, exp);
    else
        fmpz_fdiv_q_2exp(z, x, exp);
}

static inline void
_fmpz_set_si_small(fmpz_t x, slong v)
{
    fmpz_clear(x);
    *x = v;
}

slong _fmpz_sub_small_large(const fmpz_t x, const fmpz_t y);

static inline slong
_fmpz_sub_small(const fmpz_t x, const fmpz_t y)
{
    if (!COEFF_IS_MPZ(*x) && !COEFF_IS_MPZ(*y))
    {
        return (*x) - (*y);
    }
    else
    {
        return _fmpz_sub_small_large(x, y);
    }
}

static inline mp_size_t
_fmpz_size(const fmpz_t f)
{
    fmpz d = *f;

    if (!COEFF_IS_MPZ(d))
        return 1;
    else
        return mpz_size(COEFF_TO_PTR(d));
}

static inline void
fmpz_ui_mul_ui(fmpz_t r, ulong a, ulong b)
{
    if (a < (UWORD(1) << (FLINT_BITS / 2)) && b < (UWORD(1) << (FLINT_BITS / 2)))
    {
        fmpz_set_ui(r, a * b);
    }
    else
    {
        fmpz_set_ui(r, a);
        fmpz_mul_ui(r, r, b);
    }
}

static inline void
fmpz_max(fmpz_t z, const fmpz_t x, const fmpz_t y)
{
    if (fmpz_cmp(x, y) >= 0)
        fmpz_set(z, x);
    else
        fmpz_set(z, y);
}

static inline void
fmpz_min(fmpz_t z, const fmpz_t x, const fmpz_t y)
{
    if (fmpz_cmp(x, y) < 0)
        fmpz_set(z, x);
    else
        fmpz_set(z, y);
}

#define FMPZ_GET_MPN_READONLY(zsign, zn, zptr, ztmp, zv) \
    if (!COEFF_IS_MPZ(zv)) \
    { \
        (zsign) = (zv) < 0; \
        (ztmp) = FLINT_ABS(zv); \
        (zptr) = &(ztmp); \
        (zn) = 1; \
    } \
    else \
    { \
        __mpz_struct * ___zz = COEFF_TO_PTR(zv); \
        (zptr) = ___zz->_mp_d; \
        (zn) = ___zz->_mp_size; \
        (zsign) = (zn) < 0; \
        (zn) = FLINT_ABS(zn); \
    }

void fmpz_lshift_mpn(fmpz_t z, mp_srcptr d, mp_size_t dn, int sgnbit, flint_bitcnt_t shift);

static inline slong
fmpz_allocated_bytes(const fmpz_t x)
{
    if (COEFF_IS_MPZ(*x))
        return sizeof(__mpz_struct) + COEFF_TO_PTR(*x)->_mp_alloc * sizeof(mp_limb_t);
    else
        return 0;
}

#ifdef __cplusplus
}
#endif

#endif

