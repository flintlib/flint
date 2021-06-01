/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


int _fq_nmod_mpoly_compose_fq_nmod_poly_sp(fq_nmod_poly_t A, const fq_nmod_mpoly_t B,
                fq_nmod_poly_struct * const * C, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success = 1;
    flint_bitcnt_t bits = B->bits;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, shift, off;
    slong Blen = B->length;
    const mp_limb_t * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong * degrees;
    slong * offs;
    ulong * masks;
    fq_nmod_poly_struct * powers;
    fq_nmod_poly_t t, t2;
    TMP_INIT;

    FLINT_ASSERT(Blen != 0);

    TMP_START;

    degrees = TMP_ARRAY_ALLOC(nvars, slong);
    mpoly_degrees_si(degrees, Bexp, Blen, bits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_ff_poly_pow_ui_is_not_feasible(C[i]->length, degrees[i]))
        {
            success = 0;
            goto cleanup_degrees;
        }

        entries += FLINT_BIT_COUNT(degrees[i]);
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fq_nmod_poly_struct);

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
            fq_nmod_poly_init(powers + k, ctx->fqctx);
            if (j == 0)
                fq_nmod_poly_set(powers + k, C[i], ctx->fqctx);
            else
                fq_nmod_poly_mul(powers + k, powers + k - 1, powers + k - 1,
                                                                   ctx->fqctx);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fq_nmod_poly_zero(A, ctx->fqctx);
    fq_nmod_poly_init(t, ctx->fqctx);
    fq_nmod_poly_init(t2, ctx->fqctx);
    for (i = 0; i < Blen; i++)
    {
        fq_nmod_poly_fit_length(t, 1, ctx->fqctx);
        n_fq_get_fq_nmod(t->coeffs + 0, Bcoeff + d*i, ctx->fqctx);
        t->length = 1;
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fq_nmod_poly_mul(t2, t, powers + k, ctx->fqctx);
                fq_nmod_poly_swap(t, t2, ctx->fqctx);
            }
        }
        fq_nmod_poly_add(A, A, t, ctx->fqctx);
    }
    fq_nmod_poly_clear(t, ctx->fqctx);
    fq_nmod_poly_clear(t2, ctx->fqctx);

    for (k = 0; k < k_len; k++)
        fq_nmod_poly_clear(powers + k, ctx->fqctx);

cleanup_degrees:

    TMP_END;

    return success;
}

int _fq_nmod_mpoly_compose_fq_nmod_poly_mp(fq_nmod_poly_t A, const fq_nmod_mpoly_t B,
                fq_nmod_poly_struct * const * C, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success = 1;
    flint_bitcnt_t bits = B->bits;
    ulong l;
    slong i, k, N, nvars = ctx->minfo->nvars;
    slong entries, k_len, off;
    slong Blen = B->length;
    const mp_limb_t * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    flint_bitcnt_t * bitcounts;
    fq_nmod_poly_struct * powers;
    fq_nmod_poly_t t, t2;
    TMP_INIT;

    FLINT_ASSERT(Blen != 0);

    TMP_START;

    bitcounts = TMP_ARRAY_ALLOC(nvars, flint_bitcnt_t);
    degrees = TMP_ARRAY_ALLOC(nvars, fmpz);
    for (i = 0; i < nvars; i++)
        fmpz_init(degrees + i);

    mpoly_degrees_ffmpz(degrees, Bexp, Blen, bits, ctx->minfo);

    /* compute how many masks are needed */
    entries = 0;
    for (i = 0; i < nvars; i++)
    {
        if (_ff_poly_pow_fmpz_is_not_feasible(C[i]->length, degrees + i))
        {
            success = 0;
            goto cleanup_degrees;
        }

        bitcounts[i] = fmpz_bits(degrees + i);
        entries += bitcounts[i];
    }
    offs = TMP_ARRAY_ALLOC(entries, slong);
    masks = TMP_ARRAY_ALLOC(entries, ulong);
    powers = TMP_ARRAY_ALLOC(entries, fq_nmod_poly_struct);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* store bit masks for each power of two of the variables */
    k = 0;
    for (i = 0; i < nvars; i++)
    {
        off = mpoly_gen_offset_mp(i, bits, ctx->minfo);

        for (l = 0; l < bitcounts[i]; l++)
        {
            offs[k] = off + (l/FLINT_BITS);
            masks[k] = UWORD(1) << (l%FLINT_BITS);
            fq_nmod_poly_init(powers + k, ctx->fqctx);
            if (l == 0)
                fq_nmod_poly_set(powers + k, C[i], ctx->fqctx);
            else
                fq_nmod_poly_mul(powers + k, powers + k - 1, powers + k - 1,
                                                                   ctx->fqctx);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    fq_nmod_poly_zero(A, ctx->fqctx);
    fq_nmod_poly_init(t, ctx->fqctx);
    fq_nmod_poly_init(t2, ctx->fqctx);
    for (i = 0; i < Blen; i++)
    {
        fq_nmod_poly_fit_length(t, 1, ctx->fqctx);
        n_fq_get_fq_nmod(t->coeffs + 0, Bcoeff + d*i, ctx->fqctx);
        t->length = 1;
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                fq_nmod_poly_mul(t2, t, powers + k, ctx->fqctx);
                fq_nmod_poly_swap(t, t2, ctx->fqctx);
            }
        }
        fq_nmod_poly_add(A, A, t, ctx->fqctx);
    }
    fq_nmod_poly_clear(t, ctx->fqctx);
    fq_nmod_poly_clear(t2, ctx->fqctx);

    for (k = 0; k < k_len; k++)
        fq_nmod_poly_clear(powers + k, ctx->fqctx);

cleanup_degrees:

    for (i = 0; i < nvars; i++)
        fmpz_clear(degrees + i);

    TMP_END;

    return success;
}


int fq_nmod_mpoly_compose_fq_nmod_poly(fq_nmod_poly_t A,
                    const fq_nmod_mpoly_t B, fq_nmod_poly_struct * const * C,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        fq_nmod_poly_zero(A, ctx->fqctx);
        return 1;
    }

    if (B->bits <= FLINT_BITS)
    {
        return _fq_nmod_mpoly_compose_fq_nmod_poly_sp(A, B, C, ctx);
    }
    else
    {
        return _fq_nmod_mpoly_compose_fq_nmod_poly_mp(A, B, C, ctx);
    }
}
