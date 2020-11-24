/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"

void nmod_mpolyd_ctx_init(nmod_mpolyd_ctx_t dctx, slong nvars)
{
    slong i;

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }
}

void nmod_mpolyd_ctx_clear(nmod_mpolyd_ctx_t dctx)
{
    flint_free(dctx->perm);
}

void nmod_mpolyd_init(nmod_mpolyd_t poly, slong nvars)
{
    slong i;

    poly->nvars = nvars;
    poly->degb_alloc = nvars;
    poly->deg_bounds = (slong *) flint_malloc(poly->degb_alloc*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeff_alloc = WORD(16);
    poly->coeffs = (mp_limb_t *) flint_malloc(poly->coeff_alloc*sizeof(mp_limb_t));
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        poly->coeffs[i] = UWORD(0);
    }
}

void nmod_mpolyd_fit_length(nmod_mpolyd_t poly, slong len) {
    if (poly->coeff_alloc < len) {
/*flint_printf("realloc %wd -> %wd\n",poly->coeff_alloc, len);*/
        poly->coeffs = (mp_limb_t *) flint_realloc(poly->coeffs, len*sizeof(mp_limb_t));
        poly->coeff_alloc = len;
    }
}

void nmod_mpolyd_set_nvars(nmod_mpolyd_t poly, slong nvars) {

    poly->nvars = nvars;
    if (poly->degb_alloc < nvars) {
        poly->deg_bounds = (slong *) flint_realloc(poly->deg_bounds, nvars*sizeof(slong));
        poly->degb_alloc = nvars;
    }
}

void nmod_mpolyd_zero(nmod_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeffs[0] = UWORD(0);
}

void nmod_mpolyd_clear(nmod_mpolyd_t poly)
{
    flint_free(poly->deg_bounds);
    flint_free(poly->coeffs);
    poly->deg_bounds = NULL;
    poly->coeffs = NULL;
}


int nmod_mpolyd_set_degbounds(nmod_mpolyd_t A, slong * bounds)
{
    slong i;
    int success = 0;
    slong degb_prod;

    degb_prod = 1;
    for (i = 0; i < A->nvars; i++)
    {
        ulong hi;
        A->deg_bounds[i] = bounds[i];
        umul_ppmm(hi, degb_prod, degb_prod, A->deg_bounds[i]);
        if (hi != WORD(0) || degb_prod < 0)
        {
            goto done;
        }
    }

    success = 1;
    nmod_mpolyd_fit_length(A, degb_prod);

done:
    return success;
}

int nmod_mpolyd_set_degbounds_perm(nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx, slong * bounds)
{
    slong i;
    int success = 0;
    const slong * perm = dctx->perm;
    slong degb_prod;

    degb_prod = 1;
    for (i = 0; i < A->nvars; i++)
    {
        ulong hi;
        A->deg_bounds[i] = bounds[perm[i]];
        umul_ppmm(hi, degb_prod, degb_prod, A->deg_bounds[i]);
        if (hi != WORD(0) || degb_prod < 0)
        {
            goto done;
        }
    }

    success = 1;
    nmod_mpolyd_fit_length(A, degb_prod);

done:
    return success;
}


/*
    convert B to A assuming degree bounds have been set in A
*/
void nmod_mpoly_convert_to_nmod_mpolyd_degbound(nmod_mpolyd_t A,
                                const nmod_mpolyd_ctx_t dctx,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    ulong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    FLINT_ASSERT(A->nvars == nvars);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        degb_prod *= A->deg_bounds[i];
    }

    for (i = 0; i < degb_prod; i++)
    {
        A->coeffs[i] = UWORD(0);
    }

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        slong off;

        mpoly_get_monomial_ui(exps, B->exps + N*i, B->bits, ctx->minfo);
        off = 0;
        for (j = 0; j < nvars; j++)
        {
            off = exps[perm[j]] + A->deg_bounds[j]*off;
        }
        A->coeffs[off] =  B->coeffs[i];
    }

    TMP_END;
}

/*
    convert B to A - sets degree bounds in A
*/
void nmod_mpoly_convert_to_nmod_mpolyd(
                            nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx,
                          const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    slong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    nmod_mpolyd_set_nvars(A, ctx->minfo->nvars);

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (B->length == 0)
    {
        nmod_mpolyd_zero(A);
        return;
    }

    TMP_START;
    exps = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

    nmod_mpoly_degrees_si(exps, B, ctx);
    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        A->deg_bounds[i] = exps[perm[i]] + 1;
        degb_prod *= A->deg_bounds[i];
    }

    nmod_mpolyd_fit_length(A, degb_prod);
    for (i = 0; i < degb_prod; i++)
    {
        A->coeffs[i] = UWORD(0);
    }

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        slong off = 0;

        mpoly_get_monomial_ui((ulong *)exps, B->exps + N*i, B->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            off = exps[perm[j]] + A->deg_bounds[j]*off;
        }

        A->coeffs[off] = B->coeffs[i];
    }

    TMP_END;
}

/*
    Convert B to A
*/
void nmod_mpoly_convert_from_nmod_mpolyd(
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    const nmod_mpolyd_t B,
    const nmod_mpolyd_ctx_t dctx)
{
    slong off, j, k, N;
    slong bits, nvars = ctx->minfo->nvars;
    slong Alen;
    slong * perm = dctx->perm;
    slong perm_nontrivial = 0;
    ulong topmask;
    ulong * exps, * pcurexp, * pexps;
    TMP_INIT;

    FLINT_ASSERT(nvars == B->nvars);

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    /* find bits needed for the result */
    off = 1;
    for (j = 0; j < nvars; j++)
    {
        off *= B->deg_bounds[j];
        exps[perm[j]] = B->deg_bounds[j] - 1;
        perm_nontrivial |= j ^ perm[j];
    }

    FLINT_ASSERT(off <= B->coeff_alloc);

    bits = mpoly_exp_bits_required_ui(exps, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* we are going to push back terms manually */
    nmod_mpoly_fit_length_reset_bits(A, 0, bits, ctx);
    Alen = 0;

    /* find exponent vector for all variables */
    pexps = (ulong *) TMP_ALLOC(N*nvars*sizeof(ulong));
    for (k = 0; k < nvars; k++)
    {
        for (j = 0; j < nvars; j++)
            exps[perm[j]] = (j == k);
        mpoly_set_monomial_ui(pexps + k*N, exps, bits, ctx->minfo);
    }

    /* get most significant exponent in exps and its vector in ptempexp */
    off--;
    pcurexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(pcurexp, N);
    k = off;
    for (j = nvars - 1; j >= 0; j--) 
    {
        exps[j] = k % B->deg_bounds[j];
        k = k / B->deg_bounds[j];
        mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
    }

    /* scan down through the exponents */
    topmask = 0;
    for (; off >= 0; off--)
    {
        if (B->coeffs[off] != UWORD(0))
        {
            _nmod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                                   &A->exps, &A->exps_alloc, N, Alen + 1);
            A->coeffs[Alen] = B->coeffs[off];
            mpoly_monomial_set(A->exps + N*Alen, pcurexp, N);
            topmask |= (A->exps + N*Alen)[N - 1];
            Alen++;
        }

        j = nvars - 1;
        do {
            --exps[j];
            if ((slong)(exps[j]) < WORD(0))
            {
                FLINT_ASSERT(off == 0 || j > 0);
                FLINT_ASSERT(exps[j] == -UWORD(1));
                exps[j] = B->deg_bounds[j] - 1;
                mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
            } else
            {
                mpoly_monomial_sub_mp(pcurexp, pcurexp, pexps + N*j, N);
                break;
            }
        } while (--j >= 0);
    }
    _nmod_mpoly_set_length(A, Alen, ctx);

    /* sort the exponents if needed */
    if (ctx->minfo->ord != ORD_LEX || perm_nontrivial != WORD(0))
    {
        slong msb;
        mpoly_get_cmpmask(pcurexp, N, bits, ctx->minfo);
        if (topmask != WORD(0))
        {
            count_leading_zeros(msb, topmask);
            msb = (FLINT_BITS - 1)^msb;
        } else
        {
            msb = -WORD(1);
        }
        if (N == 1) {
            if (msb >= WORD(0))
            {
                _nmod_mpoly_radix_sort1(A, 0, A->length,
                                                   msb, pcurexp[0], topmask);
            }
        } else {
            _nmod_mpoly_radix_sort(A, 0, A->length,
                                        (N - 1)*FLINT_BITS + msb, N, pcurexp);
        }
    }

    TMP_END;
}


void nmod_mpolyd_print(nmod_mpolyd_t poly)
{

    int first = 0;
    slong i, j;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < poly->nvars; j++) {
        degb_prod *= poly->deg_bounds[j];
    }

    first = 1;
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (poly->coeffs[i] == 0)
            continue;

        if (!first)
            printf(" + ");

        flint_printf("%wu", poly->coeffs[i]);

        for (j = poly->nvars - 1; j >= 0; j--) 
        {
            ulong m = poly->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            flint_printf("*x%wd^%wd", j, e);
        }
        FLINT_ASSERT(k == 0);
        first = 0;
    }

    if (first)
        flint_printf("0");
}

slong nmod_mpolyd_length(const nmod_mpolyd_t A)
{
    slong i, j, degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
        degb_prod *= A->deg_bounds[j];

    for (i = degb_prod; i > 0; i--)
    {
        if (A->coeffs[i - 1] != UWORD(0))
            break;
    }

    return i;
}

