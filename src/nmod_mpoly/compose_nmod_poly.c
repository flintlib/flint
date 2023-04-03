/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


int _ff_poly_pow_fmpz_is_not_feasible(slong length, const fmpz_t e)
{
    if (length < 2)
    {
        return 0;
    }
    else
    {
        ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
        return fmpz_cmp_ui(e, limit/(ulong)(length)) >= 0;
    }
}

int _ff_poly_pow_ui_is_not_feasible(slong length, ulong e)
{
    if (length < 2)
    {
        return 0;
    }
    else
    {
        ulong limit = (ulong)(WORD_MAX)/(ulong)(2*sizeof(fmpz));
        return e >= limit/(ulong)(length);
    }
}

int _nmod_mpoly_compose_nmod_poly_sp(nmod_poly_t A, const nmod_mpoly_t B,
                      nmod_poly_struct * const * C, const nmod_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = B->bits;
    slong i, j, k, N, nvars = ctx->minfo->nvars;
    slong shift, off;
    slong entries, k_len;
    slong Blen = B->length;
    const mp_limb_t * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong * degrees;
    slong * offs;
    ulong * masks;
    nmod_poly_struct * powers;
    nmod_poly_t t, t2;
    TMP_INIT;

    FLINT_ASSERT(Blen != 0);

    TMP_START;

    degrees = (slong *) TMP_ALLOC(nvars*sizeof(slong));
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
    powers = TMP_ARRAY_ALLOC(entries, nmod_poly_struct);

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
            nmod_poly_init_mod(powers + k, ctx->mod);
            if (j == 0)
                nmod_poly_set(powers + k, C[i]);
            else
                nmod_poly_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    nmod_poly_zero(A);
    nmod_poly_init_mod(t, ctx->mod);
    nmod_poly_init_mod(t2, ctx->mod);
    for (i = 0; i < Blen; i++)
    {
        nmod_poly_zero(t);
        nmod_poly_set_coeff_ui(t, 0, Bcoeff[i]);
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                nmod_poly_mul(t2, t, powers + k);
                nmod_poly_swap(t, t2);
            }
        }
        nmod_poly_add(A, A, t);
    }
    nmod_poly_clear(t);
    nmod_poly_clear(t2);

    for (k = 0; k < k_len; k++)
        nmod_poly_clear(powers + k);

cleanup_degrees:

    TMP_END;

    return success;
}


int _nmod_mpoly_compose_nmod_poly_mp(nmod_poly_t A, const nmod_mpoly_t B,
                      nmod_poly_struct * const * C, const nmod_mpoly_ctx_t ctx)
{
    int success = 1;
    flint_bitcnt_t bits = B->bits;
    ulong l;
    slong i, k, N, nvars = ctx->minfo->nvars;
    slong off, entries, k_len;
    slong Blen = B->length;
    const mp_limb_t * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    fmpz * degrees;
    slong * offs;
    ulong * masks;
    flint_bitcnt_t * bitcounts;
    nmod_poly_struct * powers;
    nmod_poly_t t, t2;
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
    powers = TMP_ARRAY_ALLOC(entries, nmod_poly_struct);

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
            nmod_poly_init_mod(powers + k, ctx->mod);
            if (l == 0)
                nmod_poly_set(powers + k, C[i]);
            else
                nmod_poly_mul(powers + k, powers + k - 1, powers + k - 1);
            k++;
        }
    }
    k_len = k;
    FLINT_ASSERT(k_len == entries);

    /* accumulate answer */
    nmod_poly_zero(A);
    nmod_poly_init_mod(t, ctx->mod);
    nmod_poly_init_mod(t2, ctx->mod);
    for (i = 0; i < Blen; i++)
    {
        nmod_poly_zero(t);
        nmod_poly_set_coeff_ui(t, 0, Bcoeff[i]);
        for (k = 0; k < k_len; k++)
        {
            if ((Bexp[N*i + offs[k]] & masks[k]) != WORD(0))
            {
                nmod_poly_mul(t2, t, powers + k);
                nmod_poly_swap(t, t2);
            }
        }
        nmod_poly_add(A, A, t);
    }
    nmod_poly_clear(t);
    nmod_poly_clear(t2);

    for (k = 0; k < k_len; k++)
        nmod_poly_clear(powers + k);

cleanup_degrees:

    for (i = 0; i < nvars; i++)
        fmpz_clear(degrees + i);

    TMP_END;

    return success;
}


int nmod_mpoly_compose_nmod_poly(nmod_poly_t A,
                        const nmod_mpoly_t B, nmod_poly_struct * const * C,
                                                    const nmod_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        nmod_poly_zero(A);
        return 1;
    }

    if (B->bits <= FLINT_BITS)
    {
        return _nmod_mpoly_compose_nmod_poly_sp(A, B, C, ctx);
    }
    else
    {
        return _nmod_mpoly_compose_nmod_poly_mp(A, B, C, ctx);
    }
}

