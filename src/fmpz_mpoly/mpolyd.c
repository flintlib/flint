/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

#if 0
void fmpz_mpolyd_swap(fmpz_mpolyd_t A, fmpz_mpolyd_t B)
{
    fmpz_mpolyd_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mpolyd_ctx_init(fmpz_mpolyd_ctx_t dctx, slong nvars)
{
    slong i;

    dctx->nvars = nvars;
    dctx->perm = (slong *) flint_malloc(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }
}

int fmpz_mpolyd_ctx_init_version1(fmpz_mpolyd_ctx_t dctx,
                            const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success = 0;
    slong i, j, degb_prod;
    slong * Aexps, * Bexps, * deg_bounds;
    slong nvars = ctx->minfo->nvars;
    slong * perm = dctx->perm;
    TMP_INIT;

    TMP_START;
    Aexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Bexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        goto cleanup;
    fmpz_mpoly_degrees_si(Aexps, A, ctx);
    fmpz_mpoly_degrees_si(Bexps, B, ctx);

    deg_bounds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    degb_prod = 1;
    for (i = 0; i < nvars; i++)
    {
        ulong hi;
        deg_bounds[i] = FLINT_MAX(Aexps[i] + 1, Bexps[i] + 1);
        umul_ppmm(hi, degb_prod, degb_prod, deg_bounds[i]);
        if (hi != WORD(0) || degb_prod < 0)
            goto cleanup;
    }

    success = 1;
    for (i = 1; i < nvars; i++)
    {
        for (j = i; (j > 0) && (deg_bounds[j-1] < deg_bounds[j-0]); j--)
        {
            slong t1, t2;
            t1 = deg_bounds[j-1];
            t2 = deg_bounds[j-0];
            deg_bounds[j-0] = t1;
            deg_bounds[j-1] = t2;
            t1 = perm[j-1];
            t2 = perm[j-0];
            perm[j-0] = t1;
            perm[j-1] = t2;
        }
    }

cleanup:
    TMP_END;
    return success;
}

void fmpz_mpolyd_ctx_clear(fmpz_mpolyd_ctx_t dctx)
{
    flint_free(dctx->perm);
}

void fmpz_mpolyd_zero(fmpz_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    poly->coeffs[0] = UWORD(0);
}

void fmpz_mpolyd_set_fmpz(fmpz_mpolyd_t poly, fmpz_t num)
{
    slong i;

    for (i = 0; i < poly->nvars; i++)
    {
        poly->deg_bounds[i] = WORD(1);
    }
    fmpz_set(poly->coeffs + 0, num);
}

void fmpz_mpolyd_set_nvars(fmpz_mpolyd_t poly, slong nvars)
{
    poly->nvars = nvars;
    if (poly->degb_alloc < nvars) {
        poly->deg_bounds = (slong *) flint_realloc(poly->deg_bounds, nvars*sizeof(slong));
        poly->degb_alloc = nvars;
    }
}

void fmpz_mpoly_convert_to_fmpz_mpolyd(
                            fmpz_mpolyd_t poly1, const fmpz_mpolyd_ctx_t dctx,
                          const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    slong degb_prod;
    slong i, j, N;
    slong * exps;
    const slong * perm = dctx->perm;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    fmpz_mpolyd_set_nvars(poly1, ctx->minfo->nvars);

    FLINT_ASSERT(poly2->bits <= FLINT_BITS);

    if (poly2->length == 0)
    {
        fmpz_mpolyd_zero(poly1);
        return;
    }

    TMP_START;
    exps = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

    fmpz_mpoly_degrees_si(exps, poly2, ctx);
    degb_prod = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        poly1->deg_bounds[i] = exps[perm[i]] + 1;
        degb_prod *= poly1->deg_bounds[i];
    }

    fmpz_mpolyd_fit_length(poly1, degb_prod);
    for (i = 0; i < degb_prod; i++)
    {
        fmpz_zero(poly1->coeffs + i);
    }

    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);
    for (i = 0; i < poly2->length; i++)
    {
        slong off = 0;

        mpoly_get_monomial_ui((ulong *)exps, poly2->exps + N*i, poly2->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            off = exps[perm[j]] + poly1->deg_bounds[j]*off;
        }

        fmpz_set(poly1->coeffs + off, poly2->coeffs + i);
    }

    TMP_END;
}

void fmpz_mpoly_convert_from_fmpz_mpolyd(
                                  fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx,
                           const fmpz_mpolyd_t B, const fmpz_mpolyd_ctx_t dctx)
{
    slong i, j;
    slong degb_prod;
    slong * perm = dctx->perm;
    ulong * exps;
    TMP_INIT;

    FLINT_ASSERT(ctx->minfo->nvars == B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++) {
        degb_prod *= B->deg_bounds[j];
    }

    TMP_START;
    exps = (ulong *) TMP_ALLOC(B->nvars*sizeof(ulong));

    fmpz_mpoly_zero(A, ctx);
    for (i = 0; i < degb_prod; i++) {
        ulong k = i;

        if (fmpz_is_zero(B->coeffs + i))
            continue;

        for (j = B->nvars - 1; j >= 0; j--)
        {
            ulong m = B->deg_bounds[j];
            ulong e = k % m;
            k = k / m;
            exps[perm[j]] = e;
        }
        FLINT_ASSERT(k == 0);

        fmpz_mpoly_set_coeff_fmpz_ui(A, B->coeffs + i, exps, ctx);
    }

    TMP_END;
}
#endif

void fmpz_mpolyd_init(fmpz_mpolyd_t poly, slong nvars)
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
    poly->coeffs = (fmpz *) flint_malloc(poly->coeff_alloc*sizeof(fmpz));
    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fmpz_init(poly->coeffs + i);
    }
}


void fmpz_mpolyd_fit_length(fmpz_mpolyd_t poly, slong len)
{
    if (poly->coeff_alloc < len) {
        slong i;
        poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, len*sizeof(fmpz));
        for (i = poly->coeff_alloc; i < len; i++)
        {
            fmpz_init(poly->coeffs + i);
        }
        poly->coeff_alloc = len;
    }
}

void fmpz_mpolyd_clear(fmpz_mpolyd_t poly)
{
    slong i;

    for (i = 0; i < poly->coeff_alloc; i++)
    {
        fmpz_clear(poly->coeffs + i);
    }

    flint_free(poly->deg_bounds);
    flint_free(poly->coeffs);
    poly->deg_bounds = NULL;
    poly->coeffs = NULL;
}
