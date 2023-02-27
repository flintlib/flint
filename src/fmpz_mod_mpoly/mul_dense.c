/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "long_extras.h"


void _fmpz_mod_mpoly_init_dense_mock(
    fmpz_mod_poly_t D,
    const fmpz_mod_mpoly_t A,
    const slong * Adeg_bounds,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N, i, j, off, Ddeg, degb_prod;
    ulong * exps;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    degb_prod = 1;
    for (i = 0; i < nvars; i++)
        degb_prod *= Adeg_bounds[i];

    D->alloc = degb_prod;
    D->coeffs = (fmpz *) flint_calloc(degb_prod, sizeof(fmpz));

    TMP_START;
    exps = TMP_ARRAY_ALLOC(nvars, ulong);
    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);

    Ddeg = -1;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ui(exps, A->exps + N*i, A->bits, ctx->minfo);
        off = exps[0];
        for (j = 1; j < nvars; j++)
            off = exps[j] + Adeg_bounds[j]*off;

        D->coeffs[off] = A->coeffs[i]; /* shallow copy */
        Ddeg = FLINT_MAX(Ddeg, off);
    }

    D->length = Ddeg + 1;

    TMP_END;
}


static void _from_dense(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const slong * Adeg_bounds,
    fmpz_mod_poly_t D,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong off, Alen, j, k, N;
    flint_bitcnt_t bits;
    slong nvars = ctx->minfo->nvars;
    ulong topmask;
    slong * exps;
    ulong * pcurexp, * pexps;
    TMP_INIT;

    TMP_START;
    exps = TMP_ARRAY_ALLOC(nvars, slong);

    /* find bits needed for the result */
    off = 1;
    for (j = 0; j < nvars; j++)
    {
        off *= Adeg_bounds[j];
        exps[j] = Adeg_bounds[j] - 1;
    }

    bits = mpoly_exp_bits_required_ui((ulong *)exps, ctx->minfo);
    bits = FLINT_MAX(bits, Abits);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    pcurexp = TMP_ARRAY_ALLOC(N*(nvars + 1), ulong);
    pexps = pcurexp + N;

    /* we are going to push back terms manually */
    fmpz_mod_mpoly_fit_length_reset_bits(A, 0, bits, ctx);
    Alen = 0;

    /* find exponent vector for all variables */
    for (k = 0; k < nvars; k++)
        mpoly_gen_monomial_sp(pexps + k*N, k, bits, ctx->minfo);

    /* get most significant exponent in exps and its vector in ptempexp */
    off--;
    mpoly_monomial_zero(pcurexp, N);
    k = off;
    for (j = nvars - 1; j >= 0; j--) 
    {
        exps[j] = k % Adeg_bounds[j];
        k = k / Adeg_bounds[j];
        mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
    }

    /* scan down through the exponents */
    topmask = 0;
    for (; off >= 0; off--)
    {
        if (off < D->length && !fmpz_is_zero(D->coeffs + off))
        {
            _fmpz_mod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                                       &A->exps, &A->exps_alloc, N, Alen + 1);
            fmpz_swap(A->coeffs + Alen, D->coeffs + off);
            mpoly_monomial_set(A->exps + N*Alen, pcurexp, N);
            topmask |= (A->exps + N*Alen)[N - 1];
            Alen++;
        }

        j = nvars - 1;
        do {
            --exps[j];
            if (exps[j] < 0)
            {
                FLINT_ASSERT(off == 0 || j > 0);
                FLINT_ASSERT(exps[j] == -UWORD(1));
                exps[j] = Adeg_bounds[j] - 1;
                mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
            }
            else
            {
                mpoly_monomial_sub_mp(pcurexp, pcurexp, pexps + N*j, N);
                break;
            }
        } while (--j >= 0);
    }
    _fmpz_mod_mpoly_set_length(A, Alen, ctx);

    /* sort the exponents if needed */
    if (ctx->minfo->ord != ORD_LEX)
    {
        flint_bitcnt_t pos;
        mpoly_get_cmpmask(pcurexp, N, bits, ctx->minfo);
        pos = FLINT_BIT_COUNT(topmask);
        if (N == 1)
            _fmpz_mod_mpoly_radix_sort1(A->coeffs, A->exps, 0, A->length,
                                                     pos, pcurexp[0], topmask);
        else
            _fmpz_mod_mpoly_radix_sort(A->coeffs, A->exps, 0, A->length,
                                         (N - 1)*FLINT_BITS + pos, N, pcurexp);
    }

    TMP_END;
}

int _fmpz_mod_mpoly_mul_dense_maxfields(
    fmpz_mod_mpoly_t P,
    const fmpz_mod_mpoly_t A, fmpz * maxAfields,
    const fmpz_mod_mpoly_t B, fmpz * maxBfields,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i, nvars = ctx->minfo->nvars;
    fmpz_mod_poly_t Ad, Bd, Pd;
    slong prod_deg, * Abounds, * Bbounds, * Pbounds;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(nvars > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    TMP_START;

    /*
        for each variable v except for the outermost variable,
        we need to pack to degree deg_v(A) + deg_v(B)
    */
    Abounds = TMP_ARRAY_ALLOC(3*nvars, slong);
    Bbounds = Abounds + nvars;
    Pbounds = Bbounds + nvars;

    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Abounds, maxAfields, ctx->minfo);
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Bbounds, maxBfields, ctx->minfo);

    prod_deg = 1;
    for (i = 0; i < nvars; i++)
    {
        if (z_add_checked(&Abounds[i], Abounds[i], 1) ||
            z_add_checked(&Bbounds[i], Bbounds[i], 1) ||
            z_add_checked(&Pbounds[i], Abounds[i], Bbounds[i] - 1) ||
            z_mul_checked(&prod_deg, prod_deg, Pbounds[i]))
        {
            success = 0;
            goto cleanup;
        }

        if (i != 0)
        {
            /* variable of index i is not the outermost */
            Abounds[i] = Pbounds[i];
            Bbounds[i] = Pbounds[i];
        }
    }

    _fmpz_mod_mpoly_init_dense_mock(Ad, A, Abounds, ctx);
    _fmpz_mod_mpoly_init_dense_mock(Bd, B, Bbounds, ctx);
    fmpz_mod_poly_init(Pd, ctx->ffinfo);
    fmpz_mod_poly_mul(Pd, Ad, Bd, ctx->ffinfo);
    _from_dense(P, FLINT_MAX(A->bits, B->bits), Pbounds, Pd, ctx);
    fmpz_mod_poly_clear(Pd, ctx->ffinfo);
    _fmpz_mod_mpoly_clear_dense_mock(Ad);
    _fmpz_mod_mpoly_clear_dense_mock(Bd);

    success = 1;

cleanup:

    TMP_END;

    return success;
}


int fmpz_mod_mpoly_mul_dense(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                      const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length < 1 || C->length < 1)
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return 1;
    }

    if (B->bits > FLINT_BITS || C->bits > FLINT_BITS ||
        ctx->minfo->nvars < 1)
    {
        return 0;
    }

    TMP_START;

    maxBfields = TMP_ARRAY_ALLOC(2*ctx->minfo->nfields, fmpz);
    maxCfields = maxBfields + ctx->minfo->nfields;
    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    success = _fmpz_mod_mpoly_mul_dense_maxfields(A, B, maxBfields, C, maxCfields, ctx);

    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;

    return success;
}
