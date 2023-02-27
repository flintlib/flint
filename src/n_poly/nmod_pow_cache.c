/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"

/* hold positive and negative powers of b */
void nmod_pow_cache_start(
    mp_limb_t b,
    n_poly_t pos,       /* b^0, b^1, b^2, ..., b^50 */
    n_poly_t bin,       /* b^1, b^2, b^3,  b^4, b^8, b^12, ... */
    n_poly_t neg)       /* b^-0, b^-1, b^-2, ..., b^-50 */
{
    n_poly_fit_length(pos, 2);
    pos->length = 2;
    pos->coeffs[0] = 1;
    pos->coeffs[1] = b;
    bin->length = 0;
    neg->length = 0;
}

/* return a*b^e */
static mp_limb_t nmod_pow_cache_mulpow_ui_array_bin(
    mp_limb_t a,
    mp_limb_t * elimbs, slong elen,
    n_poly_t bin,
    mp_limb_t b,
    nmod_t ctx)
{
    slong ei = 0, i = 0;
    mp_limb_t e = (ei < elen) ? elimbs[ei] : 0;
    int bits_left = FLINT_BITS;

    /* complicated code needed if an odd number of bits per limb */
    FLINT_ASSERT((FLINT_BITS % 2) == 0);

    if (bin->length < 3)
    {
        n_poly_fit_length(bin, 3);
        bin->length = 3;
        bin->coeffs[0] = b;
        bin->coeffs[1] = nmod_mul(b, b, ctx);
        bin->coeffs[2] = nmod_mul(bin->coeffs[1], b, ctx);
    }

    while (ei < elen)
    {
        FLINT_ASSERT(i <= bin->length);
        if (i + 3 > bin->length)
        {
            FLINT_ASSERT(i >= 3);
            n_poly_fit_length(bin, bin->length + 3);
            bin->length += 3;
            b = nmod_mul(bin->coeffs[i - 2], bin->coeffs[i - 2], ctx);
            bin->coeffs[i + 0] = b;
            bin->coeffs[i + 1] = nmod_mul(b, b, ctx);
            bin->coeffs[i + 2] = nmod_mul(bin->coeffs[i + 1], b, ctx);
        }

        if ((e%4) != 0)
            a = nmod_mul(a, bin->coeffs[i + (e%4) - 1], ctx);

        i += 3;

        e = e/4;

        if (ei + 1 < elen)
        {
            bits_left -= 2;
            if (bits_left <= 0)
            {
                ei++;
                e = elimbs[ei];
                bits_left = FLINT_BITS;
            }
        }
        else
        {
            if (e == 0)
                break;
        }

    }

    return a;
}

/* return a*b^e */
mp_limb_t nmod_pow_cache_mulpow_ui(
    mp_limb_t a,
    ulong e,
    n_poly_t pos,
    n_poly_t bin,
    n_poly_t neg,
    nmod_t ctx)
{
    slong i;
    mp_limb_t b;

    FLINT_ASSERT(pos->length >= 2);

    b = pos->coeffs[1];
    if (b <= 1)
        return (b == 1 || e == 0) ? a : 0;

    if (e < 50)
    {
        n_poly_fit_length(pos, e + 1);
        i = pos->length;
        while (i <= e)
        {
            pos->coeffs[i] = nmod_mul(b, pos->coeffs[i - 1], ctx);
            pos->length = ++i;
        }

        return nmod_mul(a, pos->coeffs[e], ctx);
    }

    return nmod_pow_cache_mulpow_ui_array_bin(a, &e, 1, bin, b, ctx);
}

/* return a*b^-e, assume ctx.n is prime */
mp_limb_t nmod_pow_cache_mulpow_neg_ui(
    mp_limb_t a,
    ulong e,
    n_poly_t pos,
    n_poly_t bin,
    n_poly_t neg,
    nmod_t ctx)
{
    slong i;
    mp_limb_t b;

    FLINT_ASSERT(pos->length >= 2);

    b = pos->coeffs[1];
    if (b <= 1)
        return (b == 1 || e == 0) ? a : 0;

    if (e < 50)
    {
        if (neg->length < 2)
        {
            n_poly_fit_length(neg, 2);
            neg->length = 2;
            neg->coeffs[0] = 1;
            neg->coeffs[1] = nmod_inv(b, ctx);            
        }

        n_poly_fit_length(neg, e + 1);

        i = neg->length;
        while (i <= e)
        {
            neg->coeffs[i] = nmod_mul(neg->coeffs[1],
                                             neg->coeffs[i - 1], ctx);
            neg->length = ++i;
        }

        return nmod_mul(a, neg->coeffs[e], ctx);
    }

    if (e >= ctx.n)
        e = e % (ctx.n - 1);
    e = ctx.n - 1 - e;

    return nmod_pow_cache_mulpow_ui(a, e, pos, bin, neg, ctx);
}

/* return a*b^-e */
mp_limb_t nmod_pow_cache_mulpow_fmpz(
    mp_limb_t a,
    const fmpz_t e,
    n_poly_t pos,
    n_poly_t bin,
    n_poly_t neg,
    nmod_t ctx)
{
    mp_limb_t b = pos->coeffs[1];

    FLINT_ASSERT(pos->length >= 2);

    if (b <= 1)
        return (b == 1 || fmpz_is_zero(e)) ? a : 0;

    if (!COEFF_IS_MPZ(*e))
    {
        if (*e >= 0)
            return nmod_pow_cache_mulpow_ui(a, *e, pos, bin, neg, ctx);
        else
            return nmod_pow_cache_mulpow_neg_ui(a, -*e, pos, bin, neg, ctx);
    }
    else
    {
        if (COEFF_TO_PTR(*e)->_mp_size >= 0)
            return nmod_pow_cache_mulpow_ui_array_bin(a, COEFF_TO_PTR(*e)->_mp_d,
                                      COEFF_TO_PTR(*e)->_mp_size, bin, b, ctx);
        else
            return nmod_pow_cache_mulpow_ui(a, fmpz_fdiv_ui(e, ctx.n - 1),
                                                           pos, bin, neg, ctx);
    }
}
