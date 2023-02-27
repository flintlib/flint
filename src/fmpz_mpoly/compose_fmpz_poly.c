/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

static int _fmpz_poly_pow_fmpz_is_not_feasible(const fmpz_poly_t b, const fmpz_t e)
{
    if (b->length < 2)
    {
        if (b->length == 1)
        {
            return _fmpz_pow_fmpz_is_not_feasible(fmpz_bits(b->coeffs + 0), e);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
        return fmpz_cmp_ui(e, limit/(ulong)(b->length)) >= 0;
    }
}

static int _fmpz_poly_pow_ui_is_not_feasible(const fmpz_poly_t b, ulong e)
{
    if (b->length < 2)
    {
        if (b->length == 1)
        {
            return _fmpz_pow_ui_is_not_feasible(fmpz_bits(b->coeffs + 0), e);
        }
        else
        {
            return 0;
        }
    }
    else
    {
        ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
        return e >= limit/(ulong)(b->length);
    }
}

int _fmpz_mpoly_compose_fmpz_poly_sp(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = B->bits;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, shift, off;
    slong Blen = B->length;
    const fmpz * Bcoeff = B->coeffs;
    const ulong * Bexp = B->exps;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    fmpz_poly_struct * powers;
    fmpz_poly_t t, t2;
    TMP_INIT;

    FLINT_ASSERT(Blen != 0);

    TMP_START;

    degrees = TMP_ARRAY_ALLOC(nvars, slong);
    mpoly_degrees_si(degrees, Bexp, Blen, bits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_fmpz_poly_pow_ui_is_not_feasible(C[i], degrees[i]))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += FLINT_BIT_COUNT(degrees[i]);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fmpz_poly_struct);

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
            fmpz_poly_init(powers + k);
            if (j == 0)
                fmpz_poly_set(powers + k, C[i]);
            else
                fmpz_poly_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fmpz_poly_zero(A);
    fmpz_poly_init(t);
    fmpz_poly_init(t2);
    for (i = 0; i < Blen; i++)
    {
        fmpz_poly_set_fmpz(t, Bcoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fmpz_poly_mul(t2, t, powers + k);
                fmpz_poly_swap(t, t2);
            }
        }
        fmpz_poly_add(A, A, t);
    }
    fmpz_poly_clear(t);
    fmpz_poly_clear(t2);

    for (k = 0; k < k_len; k++)
        fmpz_poly_clear(powers + k);

cleanup_degrees:

    TMP_END;

    return success;
}


int _fmpz_mpoly_compose_fmpz_poly_mp(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = B->bits;
    ulong l;
    slong i, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, off;
    slong Blen = B->length;
    fmpz * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    fmpz_poly_struct * powers;
    fmpz_poly_t t, t2;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);

    TMP_START;

    degrees = TMP_ARRAY_ALLOC(nvars, fmpz);
    for (i = 0; i < nvars; i++)
        fmpz_init(degrees + i);

    mpoly_degrees_ffmpz(degrees, Bexp, Blen, bits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_fmpz_poly_pow_fmpz_is_not_feasible(C[i], degrees + i))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += fmpz_bits(degrees + i);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fmpz_poly_struct);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        flint_bitcnt_t varibits = fmpz_bits(degrees + i);

        off = mpoly_gen_offset_mp(i, bits, ctx->minfo);

        for (l = 0; l < varibits; l++)
        {
            offs[k] = off + (l / FLINT_BITS);
            masks[k] = UWORD(1) << (l % FLINT_BITS);
            fmpz_poly_init(powers + k);
            if (l == 0)
                fmpz_poly_set(powers + k, C[i]);
            else
                fmpz_poly_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fmpz_poly_zero(A);
    fmpz_poly_init(t);
    fmpz_poly_init(t2);
    for (i = 0; i < Blen; i++)
    {
        fmpz_poly_set_fmpz(t, Bcoeff + i);
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fmpz_poly_mul(t2, t, powers + k);
                fmpz_poly_swap(t, t2);
            }
        }
        fmpz_poly_add(A, A, t);
    }
    fmpz_poly_clear(t);
    fmpz_poly_clear(t2);

    for (k = 0; k < k_len; k++)
        fmpz_poly_clear(powers + k);

cleanup_degrees:

    for (i = 0; i < nvars; i++)
        fmpz_clear(degrees + i);

    TMP_END;

    return success;
}


int fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B,
                      fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fmpz_poly_zero(A);
        return 1;
    }

    if (B->bits <= FLINT_BITS)
    {
        return _fmpz_mpoly_compose_fmpz_poly_sp(A, B, C, ctx);
    }
    else
    {
        return _fmpz_mpoly_compose_fmpz_poly_mp(A, B, C, ctx);
    }
}

