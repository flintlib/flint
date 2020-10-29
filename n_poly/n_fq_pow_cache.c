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
void n_fq_pow_cache_start_n_fq(
    const mp_limb_t * b,
    n_poly_t pos,       /* b^0, b^1, b^2, ..., b^50 */
    n_poly_t bin,       /* b^1, b^2, b^3,  b^4, b^8, b^12, ... */
    n_poly_t neg,       /* b^-0, b^-1, b^-2, ..., b^-50 */
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    n_poly_fit_length(pos, d*2);
    pos->length = 2;
    _n_fq_one(pos->coeffs + d*0, d);
    _n_fq_set(pos->coeffs + d*1, b, d);
    bin->length = 0;
    neg->length = 0;
}

void n_fq_pow_cache_start_fq_nmod(
    const fq_nmod_t b,
    n_poly_t pos,
    n_poly_t bin,
    n_poly_t neg,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    n_poly_fit_length(pos, d*2);
    pos->length = 2;
    _n_fq_one(pos->coeffs + d*0, d);
    n_fq_set_fq_nmod(pos->coeffs + d*1, b, ctx);
    bin->length = 0;
    neg->length = 0;
}

/* r = a*b^e */
static void n_fq_pow_cache_mulpow_ui_array_bin(
    mp_limb_t * r,
    const mp_limb_t * a,
    mp_limb_t * elimbs, slong elen,
    n_poly_t bin,
    const mp_limb_t * b,
    const fq_nmod_ctx_t ctx,
    mp_limb_t * tmp)    /* size d*N_FQ_MUL_ITCH */
{
    slong d = fq_nmod_ctx_degree(ctx);
    const mp_limb_t * s = a; /* source */
    slong ei = 0, i = 0;
    mp_limb_t e = (ei < elen) ? elimbs[ei] : 0;
    int bits_left = FLINT_BITS;

    /* complicated code needed if an odd number of bits per limb */
    FLINT_ASSERT((FLINT_BITS % 2) == 0);

    if (bin->length < 3)
    {
        n_poly_fit_length(bin, 3*d);
        bin->length = 3;
        _n_fq_set(bin->coeffs + d*0, b, d);
        _n_fq_mul(bin->coeffs + d*1, bin->coeffs + d*0,
                                     bin->coeffs + d*0, ctx, tmp);
        _n_fq_mul(bin->coeffs + d*2, bin->coeffs + d*1,
                                     bin->coeffs + d*0, ctx, tmp);
    }

    while (ei < elen)
    {
        FLINT_ASSERT(i <= bin->length);

        if (i + 3 > bin->length)
        {
            FLINT_ASSERT(i >= 3);

            n_poly_fit_length(bin, d*(bin->length + 3));
            bin->length += 3;

            _n_fq_mul(bin->coeffs + d*(i + 0), bin->coeffs + d*(i - 2),
                                            bin->coeffs + d*(i - 2), ctx, tmp);

            _n_fq_mul(bin->coeffs + d*(i + 1), bin->coeffs + d*(i + 0),
                                            bin->coeffs + d*(i + 0), ctx, tmp);

            _n_fq_mul(bin->coeffs + d*(i + 2), bin->coeffs + d*(i + 1),
                                            bin->coeffs + d*(i + 0), ctx, tmp);
        }

        if ((e%4) != 0)
        {
            _n_fq_mul(r, s, bin->coeffs + d*(i + (e%4) - 1), ctx, tmp);
            s = r;
        }

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

    if (s != r)
        _n_fq_set(r, s, d);
}

/* r =  a*b^e */
void n_fq_pow_cache_mulpow_ui(
    mp_limb_t * r,
    const mp_limb_t * a,
    ulong e,
    n_poly_t pos,
    n_poly_t bin,
    n_poly_t neg,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i = pos->length;
    int a_in_fp = _n_fq_is_ui(a, d);

    FLINT_ASSERT(i >= 2);

    if (a[0] == 0 && a_in_fp)
    {
        _n_fq_zero(r, d);
        return;
    }

    if (e < 50)
    {
        n_poly_fit_length(pos, d*(FLINT_MAX(e + 1, i) + N_FQ_MUL_ITCH));
        while (i <= e)
        {
            FLINT_ASSERT(d*(i + 1 + N_FQ_MUL_ITCH) <= pos->alloc);
            _n_fq_mul(pos->coeffs + d*i, pos->coeffs + d*1,
                         pos->coeffs + d*(i - 1), ctx, pos->coeffs + d*(i + 1));
            i++;
            pos->length = i;
        }

        if (a_in_fp)
            _nmod_vec_scalar_mul_nmod(r, pos->coeffs + d*e, d, a[0], ctx->mod);
        else
            _n_fq_mul(r, a, pos->coeffs + d*e, ctx, pos->coeffs + d*i);
        return;
    }

    if (_n_fq_is_zero(pos->coeffs + d*1, d))
    {
        _n_fq_zero(r, d);
        return;
    }

    n_poly_fit_length(pos, d*(i + N_FQ_MUL_ITCH));
    n_fq_pow_cache_mulpow_ui_array_bin(r, a, &e, 1, bin, pos->coeffs + d*1,
                                                       ctx, pos->coeffs + d*i);
}

/* r =  a*b^-e */
void n_fq_pow_cache_mulpow_neg_ui(
    mp_limb_t * r,
    const mp_limb_t * a,
    ulong e,
    n_poly_t pos,
    n_poly_t bin,
    n_poly_t neg,
    const fq_nmod_ctx_t ctx)
{
    slong i, d = fq_nmod_ctx_degree(ctx);
    mp_limb_t * tmp;
    fmpz_t f;

    FLINT_ASSERT(pos->length >= 2);

    if (_n_fq_is_zero(pos->coeffs + d*1, d))
    {
        if (e > 0)
            _n_fq_zero(r, d);
        else
            _n_fq_set(r, a, d);
        return;
    }

    if (e < 50)
    {
        n_poly_fit_length(pos, d*(pos->length
                                       + FLINT_MAX(N_FQ_MUL_ITCH, N_FQ_INV_ITCH)));
        tmp = pos->coeffs + d*pos->length;

        if (neg->length < 2)
        {
            n_poly_fit_length(neg, 2*d);
            neg->length = 2;
            _n_fq_one(neg->coeffs + d*0, d);
            _n_fq_inv(neg->coeffs + d*1, pos->coeffs + d*1, ctx, tmp);
        }

        i = neg->length;

        n_poly_fit_length(neg, d*(e + 1));
        while (i <= e)
        {
            _n_fq_mul(neg->coeffs + d*i, neg->coeffs + d*1,
                                     neg->coeffs + d*(i - 1), ctx, tmp);
            i++;
            neg->length = i;
        }

        _n_fq_mul(r, a, neg->coeffs + d*e, ctx, tmp);
        return;
    }

    fmpz_init(f);
    fmpz_neg_ui(f, e);
    n_fq_pow_cache_mulpow_fmpz(r, a, f, pos, bin, neg, ctx);
    fmpz_clear(f);
}

/* r = a*b^-e */
void n_fq_pow_cache_mulpow_fmpz(
    mp_limb_t * r,
    const mp_limb_t * a,
    const fmpz_t e,
    n_poly_t pos,
    n_poly_t bin,
    n_poly_t neg,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    fmpz_t t;

    FLINT_ASSERT(pos->length >= 2);

    if (!COEFF_IS_MPZ(*e) && *e >= 0)
    {
        n_fq_pow_cache_mulpow_ui(r, a, *e, pos, bin, neg, ctx);
        return;
    }

    if (_n_fq_is_zero(pos->coeffs + d*1, d))
    {
        if (!fmpz_is_zero(e))
            _n_fq_zero(r, d);
        else
            _n_fq_set(r, a, d);
        return;
    }

    fmpz_init(t);
    fq_nmod_ctx_order(t, ctx);
    fmpz_sub_ui(t, t, 1);
    fmpz_mod(t, e, t);

    n_poly_fit_length(pos, d*(pos->length + N_FQ_MUL_ITCH));

    if (COEFF_IS_MPZ(*t))
        n_fq_pow_cache_mulpow_ui_array_bin(r, a,
                     COEFF_TO_PTR(*t)->_mp_d, COEFF_TO_PTR(*t)->_mp_size, bin,
                          pos->coeffs + d*1, ctx, pos->coeffs + d*pos->length);
    else
        n_fq_pow_cache_mulpow_ui(r, a, *t, pos, bin, neg, ctx);

    fmpz_clear(t);
}
