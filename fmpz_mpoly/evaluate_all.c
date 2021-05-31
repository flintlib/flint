/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/* given the exponent and the bit count of the base, can we expect b^e to fail */
int _fmpz_pow_fmpz_is_not_feasible(flint_bitcnt_t bbits, const fmpz_t e)
{
    ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
    FLINT_ASSERT(fmpz_sgn(e) >= 0);
    return bbits > 1 && fmpz_cmp_ui(e, limit/bbits) >= 0;
}

int _fmpz_pow_ui_is_not_feasible(flint_bitcnt_t bbits, ulong e)
{
    ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
    return bbits > 1 && e >= limit/bbits;
}

int _fmpz_mpoly_evaluate_all_fmpz_sp(fmpz_t ev, const fmpz_mpoly_t A,
                                fmpz * const * val, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = A->bits;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, shift, off;
    slong Alen = A->length;
    const fmpz * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpz * powers;
    fmpz_t t;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    TMP_START;

    degrees = TMP_ARRAY_ALLOC(nvars, slong);
    mpoly_degrees_si(degrees, Aexp, Alen, bits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_fmpz_pow_ui_is_not_feasible(fmpz_bits(val[i]), degrees[i]))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += FLINT_BIT_COUNT(degrees[i]);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fmpz);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        flint_bitcnt_t varibits = FLINT_BIT_COUNT(degrees[i]);

        mpoly_gen_offset_shift_sp(&off, &shift, i, bits, ctx->minfo);
        for (j = 0; j < varibits; j++)
        {
            offs[k] = off;
            masks[k] = UWORD(1) << (shift + j);
            fmpz_init(powers + k);
            if (j == 0)
                fmpz_set(powers + k, val[i]);
            else
                fmpz_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fmpz_zero(ev);
    fmpz_init(t);
    for (i = 0; i < Alen; i++)
    {
        fmpz_set(t, Acoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Aexp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpz_mul(t, t, powers + k);
        }
        fmpz_add(ev, ev, t);
    }
    fmpz_clear(t);

    for (k = 0; k < k_len; k++)
        fmpz_clear(powers + k);

cleanup_degrees:

    TMP_END;

    return success;
}

int _fmpz_mpoly_evaluate_all_fmpz_mp(fmpz_t ev, const fmpz_mpoly_t A,
                               fmpz * const * vals, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t Abits = A->bits;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, off;
    slong Alen = A->length;
    const fmpz * Acoeff = A->coeffs;
    const ulong * Aexp = A->exps;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fmpz * powers;
    fmpz_t t;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    TMP_START;

    degrees = _fmpz_vec_init(nvars);
    mpoly_degrees_ffmpz(degrees, Aexp, Alen, Abits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_fmpz_pow_fmpz_is_not_feasible(fmpz_bits(vals[i]), degrees + i))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += fmpz_bits(degrees + i);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fmpz);

    N = mpoly_words_per_exp(Abits, ctx->minfo);

    /* store bit masks for each power of two of the variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        flint_bitcnt_t varibits = fmpz_bits(degrees + i);

        off = mpoly_gen_offset_mp(i, Abits, ctx->minfo);
        for (j = 0; j < varibits; j++)
        {
            offs[k] = off + (j / FLINT_BITS);
            masks[k] = UWORD(1) << (j % FLINT_BITS);
            fmpz_init(powers + k);
            if (j == 0)
                fmpz_set(powers + k, vals[i]);
            else
                fmpz_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fmpz_zero(ev);
    fmpz_init(t);
    for (i = 0; i < Alen; i++)
    {
        fmpz_set(t, Acoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Aexp[N*i + offs[k]] & masks[k]) != WORD(0))
                fmpz_mul(t, t, powers + k);
        }
        fmpz_add(ev, ev, t);
    }
    fmpz_clear(t);

    for (k = 0; k < k_len; k++)
        fmpz_clear(powers + k);

cleanup_degrees:

    _fmpz_vec_clear(degrees, nvars);

    TMP_END;

    return success;
}


/* evaluate A(xbar) at xbar = vals */
int fmpz_mpoly_evaluate_all_fmpz(fmpz_t ev, const fmpz_mpoly_t A,
                               fmpz * const * vals, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx))
    {
        fmpz_zero(ev);
        return 1;
    }

    if (A->bits <= FLINT_BITS)
    {
        return _fmpz_mpoly_evaluate_all_fmpz_sp(ev, A, vals, ctx);
    }
    else
    {
        return _fmpz_mpoly_evaluate_all_fmpz_mp(ev, A, vals, ctx);
    }
}
