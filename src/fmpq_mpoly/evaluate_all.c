/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

static int _fmpz_mpoly_evaluate_all_tree_fmpq_sp(fmpq_t ev, const fmpz_mpoly_t poly,
                               fmpq * const * vals, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = poly->bits;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, shift, off;
    slong p_len = poly->length;
    const fmpz * p_coeff = poly->coeffs;
    const ulong * p_exp = poly->exps;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpq * powers;
    fmpq_t t;
    TMP_INIT;

    FLINT_ASSERT(p_len > 0);

    TMP_START;

    degrees = TMP_ARRAY_ALLOC(nvars, slong);
    fmpz_mpoly_degrees_si(degrees, poly, ctx);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_fmpz_pow_ui_is_not_feasible(fmpq_height_bits(vals[i]), degrees[i]))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += FLINT_BIT_COUNT(degrees[i]);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fmpq);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        flint_bitcnt_t varibits = FLINT_BIT_COUNT(degrees[i]);

        mpoly_gen_offset_shift_sp(&off, &shift, i, bits, ctx->minfo);
        for (j = 0; j < varibits; j++)
        {
            offs[k] = off;
            masks[k] = UWORD(1) << (shift + j);
            fmpq_init(powers + k);
            if (j == 0)
                fmpq_set(powers + k, vals[i]);
            else
                fmpq_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fmpq_zero(ev);
    fmpq_init(t);
    for (i = 0; i < p_len; i++)
    {
        fmpz_set(fmpq_numref(t), p_coeff + i);
        fmpz_one(fmpq_denref(t));
        for (k = 0; k < k_len; k++)
        {
            if ((p_exp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpq_mul(t, t, powers + k);
        }
        fmpq_add(ev, ev, t);
    }
    fmpq_clear(t);

    for (k = 0; k < k_len; k++)
        fmpq_clear(powers + k);

cleanup_degrees:

    TMP_END;

    return success;
}

static int _fmpz_mpoly_evaluate_all_fmpq_mp(fmpq_t ev, const fmpz_mpoly_t poly,
                               fmpq * const * vals, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = poly->bits;
    slong i, j, k, N,  nvars = ctx->minfo->nvars;
    slong entries, k_len, off;
    slong p_len = poly->length;
    const fmpz * p_coeff = poly->coeffs;
    ulong * p_exp = poly->exps;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpq * powers;
    fmpq_t t;
    TMP_INIT;

    FLINT_ASSERT(p_len > 0);

    TMP_START;

    degrees = _fmpz_vec_init(nvars);
    mpoly_degrees_ffmpz(degrees, p_exp, p_len, bits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_fmpz_pow_fmpz_is_not_feasible(fmpq_height_bits(vals[i]), degrees + i))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += fmpz_bits(degrees + i);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fmpq);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the non-main variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        flint_bitcnt_t varibits = fmpz_bits(degrees + i);

        off = mpoly_gen_offset_mp(i, bits, ctx->minfo);
        for (j = 0; j < varibits; j++)
        {
            offs[k] = off + (j / FLINT_BITS);
            masks[k] = UWORD(1) << (j % FLINT_BITS);
            fmpq_init(powers + k);
            if (j == 0)
                fmpq_set(powers + k, vals[i]);
            else
                fmpq_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fmpq_zero(ev);
    fmpq_init(t);
    for (i = 0; i < p_len; i++)
    {
        fmpz_set(fmpq_numref(t), p_coeff + i);
        fmpz_one(fmpq_denref(t));
        for (k = 0; k < k_len; k++)
        {
            if ((p_exp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpq_mul(t, t, powers + k);
        }
        fmpq_add(ev, ev, t);
    }
    fmpq_clear(t);

    for (k = 0; k < k_len; k++)
        fmpq_clear(powers + k);

cleanup_degrees:

    _fmpz_vec_clear(degrees, nvars);

    TMP_END;

    return success;
}


/* evaluate A(xbar) at xbar = vals */
int fmpq_mpoly_evaluate_all_fmpq(fmpq_t ev, const fmpq_mpoly_t A,
                               fmpq * const * vals, const fmpq_mpoly_ctx_t ctx)
{
    int success;
    fmpq_t t;

    if (fmpq_mpoly_is_zero(A, ctx))
    {
        fmpq_zero(ev);
        return 1;
    }

    fmpq_init(t);

    if (A->zpoly->bits <= FLINT_BITS)
    {
        success = _fmpz_mpoly_evaluate_all_tree_fmpq_sp(t, A->zpoly, vals, ctx->zctx);
    }
    else
    {
        success = _fmpz_mpoly_evaluate_all_fmpq_mp(t, A->zpoly, vals, ctx->zctx);
    }

    if (success)
        fmpq_mul(ev, t, A->content);

    fmpq_clear(t);
    return success;
}
