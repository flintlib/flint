/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


int fq_zech_mpoly_pfrac_init(
    fq_zech_mpoly_pfrac_t I,
    flint_bitcnt_t bits,
    slong r,
    slong w,
    const fq_zech_mpoly_struct * betas,
    const fq_zech_struct * alpha,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong success = 1;
    slong i, j, k;
    fq_zech_poly_t p;
    fq_zech_poly_t G, S, pq;

    FLINT_ASSERT(bits <= FLINT_BITS);

    I->bits = bits;
    I->r = r;
    I->w = w;

    I->dbetas = (fq_zech_poly_struct *) flint_malloc(
                                                r*sizeof(fq_zech_poly_struct));

    I->inv_prod_dbetas = (fq_zech_poly_struct *) flint_malloc(
                                               r* sizeof(fq_zech_poly_struct));

    I->prod_mbetas = (fq_zech_mpoly_struct *) flint_malloc(
                                       (w + 1)*r*sizeof(fq_zech_mpoly_struct));

    I->prod_mbetas_coeffs = (fq_zech_mpolyv_struct *) flint_malloc(
                                      (w + 1)*r*sizeof(fq_zech_mpolyv_struct));

    I->mbetas = (fq_zech_mpoly_struct *) flint_malloc(
                                       (w + 1)*r*sizeof(fq_zech_mpoly_struct));

    I->deltas = (fq_zech_mpoly_struct *) flint_malloc(
                                       (w + 1)*r*sizeof(fq_zech_mpoly_struct));

    I->xalpha = (fq_zech_mpoly_struct *) flint_malloc(
                                         (w + 1)*sizeof(fq_zech_mpoly_struct));

    I->q = (fq_zech_mpoly_struct *) flint_malloc(
                                         (w + 1)*sizeof(fq_zech_mpoly_struct));

    I->qt = (fq_zech_mpoly_struct *) flint_malloc(
                                         (w + 1)*sizeof(fq_zech_mpoly_struct));

    I->newt = (fq_zech_mpoly_struct *) flint_malloc(
                                         (w + 1)*sizeof(fq_zech_mpoly_struct));

    I->delta_coeffs = (fq_zech_mpolyv_struct *) flint_malloc(
                                      (w + 1)*r*sizeof(fq_zech_mpolyv_struct));

    for (i = 0; i <= w; i++)
    {
        fq_zech_mpoly_init(I->xalpha + i, ctx);
        fq_zech_mpoly_init(I->q + i, ctx);
        fq_zech_mpoly_init(I->qt + i, ctx);
        fq_zech_mpoly_init(I->newt + i, ctx);
        for (j = 0; j < r; j++)
        {
            fq_zech_mpoly_init(I->deltas + i*r + j, ctx);
            fq_zech_mpolyv_init(I->delta_coeffs + i*r + j, ctx);
        }

        if (i < 1)
            continue;

        fq_zech_mpoly_gen(I->xalpha + i, i, ctx);
        fq_zech_mpoly_sub_fq_zech(I->xalpha + i, I->xalpha + i, alpha + i - 1, ctx);
        fq_zech_mpoly_repack_bits_inplace(I->xalpha + i, I->bits, ctx);
    }

    fq_zech_poly_init(p, ctx->fqctx);
    fq_zech_poly_init(G, ctx->fqctx);
    fq_zech_poly_init(S, ctx->fqctx);
    fq_zech_poly_init(pq, ctx->fqctx);

    /* set betas */
    i = w;
    for (j = 0; j < r; j++)
    {
        fq_zech_mpoly_init(I->mbetas + i*r + j, ctx);
        fq_zech_mpoly_set(I->mbetas + i*r + j, betas + j, ctx);
    }
    for (i--; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fq_zech_mpoly_init(I->mbetas + i*r + j, ctx);
            fq_zech_mpoly_evaluate_one_fq_zech(I->mbetas + i*r + j,
                             I->mbetas + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    for (j = 0; j < r; j++)
    {
        fq_zech_poly_init(I->dbetas + j, ctx->fqctx);
        if (!fq_zech_mpoly_get_fq_zech_poly(I->dbetas + j,
                                                  I->mbetas + 0*r + j, 0, ctx))
        {
            success = 0;
        }
    }

    /* set product of betas */
    for (i = w; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fq_zech_mpoly_init(I->prod_mbetas + i*r + j, ctx);
            fq_zech_mpoly_one(I->prod_mbetas + i*r + j, ctx);
            for (k = 0; k < r; k++)
            {
                if (k == j)
                    continue;
                fq_zech_mpoly_mul(I->prod_mbetas + i*r + j,
                           I->prod_mbetas + i*r + j, I->mbetas + i*r + k, ctx);
            }
            fq_zech_mpolyv_init(I->prod_mbetas_coeffs + i*r + j, ctx);
            if (i > 0)
            {
                fq_zech_mpoly_to_mpolyv(I->prod_mbetas_coeffs + i*r + j,
                                 I->prod_mbetas + i*r + j, I->xalpha + i, ctx);
            }
        }        
    }

    for (j = 0; j < r; j++)
        fq_zech_poly_init(I->inv_prod_dbetas + j, ctx->fqctx);

    for (j = 0; success && j < r; j++)
    {
        if (fq_zech_poly_degree(I->dbetas + j, ctx->fqctx) !=
            fq_zech_mpoly_degree_si(betas + j, 0, ctx))
        {
            success = 0;
        }
    }

    for (j = 0; success && j < r; j++)
    {
        fq_zech_poly_one(pq, ctx->fqctx);
        for (k = 0; k < r; k++)
        {
            if (k == j)
                continue;
            fq_zech_poly_mul(pq, pq, I->dbetas + k, ctx->fqctx);
        }
        fq_zech_poly_xgcd(G, S, I->inv_prod_dbetas + j,
                                                I->dbetas + j, pq, ctx->fqctx);
        if (!fq_zech_poly_is_one(G, ctx->fqctx))
        {
            success = 0;
        }
    }

    fq_zech_poly_clear(p, ctx->fqctx);
    fq_zech_poly_clear(G, ctx->fqctx);
    fq_zech_poly_clear(S, ctx->fqctx);
    fq_zech_poly_clear(pq, ctx->fqctx);

    I->dbetas_mvar = (fq_zech_mpoly_struct *) flint_malloc(
                                               r*sizeof(fq_zech_mpoly_struct));
    I->inv_prod_dbetas_mvar = (fq_zech_mpoly_struct *) flint_malloc(
                                               r*sizeof(fq_zech_mpoly_struct));
    for (j = 0; j < r; j++)
    {
        fq_zech_mpoly_init(I->dbetas_mvar + j, ctx);
        fq_zech_mpoly_init(I->inv_prod_dbetas_mvar + j, ctx);

        _fq_zech_mpoly_set_fq_zech_poly(I->dbetas_mvar + j, I->bits,
                             I->dbetas[j].coeffs, I->dbetas[j].length, 0, ctx);

        _fq_zech_mpoly_set_fq_zech_poly(I->inv_prod_dbetas_mvar + j, I->bits,
           I->inv_prod_dbetas[j].coeffs, I->inv_prod_dbetas[j].length, 0, ctx);
    }

    fq_zech_mpoly_init(I->T, ctx);
    fq_zech_mpoly_init(I->Q, ctx);
    fq_zech_mpoly_init(I->R, ctx);

    FLINT_ASSERT(success == 1);

    return success;
}


void fq_zech_mpoly_pfrac_clear(fq_zech_mpoly_pfrac_t I, const fq_zech_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i <= I->w; i++)
    {
        fq_zech_mpoly_clear(I->xalpha + i, ctx);
        fq_zech_mpoly_clear(I->q + i, ctx);
        fq_zech_mpoly_clear(I->qt + i, ctx);
        fq_zech_mpoly_clear(I->newt + i, ctx);
        for (j = 0; j < I->r; j++)
            fq_zech_mpolyv_clear(I->delta_coeffs + i*I->r + j, ctx);
    }

    flint_free(I->xalpha);
    flint_free(I->q);
    flint_free(I->qt);
    flint_free(I->newt);
    flint_free(I->delta_coeffs);

    for (j = 0; j < I->r; j++)
    {
        fq_zech_poly_clear(I->inv_prod_dbetas + j, ctx->fqctx);
        fq_zech_poly_clear(I->dbetas + j, ctx->fqctx);
        for (i = 0; i <= I->w; i++)
        {
            fq_zech_mpolyv_clear(I->prod_mbetas_coeffs + i*I->r + j, ctx);
            fq_zech_mpoly_clear(I->prod_mbetas + i*I->r + j, ctx);
            fq_zech_mpoly_clear(I->mbetas + i*I->r + j, ctx);
            fq_zech_mpoly_clear(I->deltas + i*I->r + j, ctx);
        }
    }
    flint_free(I->inv_prod_dbetas);
    flint_free(I->dbetas);
    flint_free(I->prod_mbetas);
    flint_free(I->prod_mbetas_coeffs);
    flint_free(I->mbetas);
    flint_free(I->deltas);

    for (j = 0; j < I->r; j++)
    {
        fq_zech_mpoly_clear(I->dbetas_mvar + j, ctx);
        fq_zech_mpoly_clear(I->inv_prod_dbetas_mvar + j, ctx);
    }
    flint_free(I->dbetas_mvar);
    flint_free(I->inv_prod_dbetas_mvar);

    fq_zech_mpoly_clear(I->T, ctx);
    fq_zech_mpoly_clear(I->Q, ctx);
    fq_zech_mpoly_clear(I->R, ctx);
}


int fq_zech_mpoly_pfrac(
    slong l,
    fq_zech_mpoly_t t,
    const slong * degs,
    fq_zech_mpoly_pfrac_t I,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, j, k;
    int success;
    fq_zech_mpoly_struct * deltas = I->deltas + l*I->r;
    fq_zech_mpoly_struct * newdeltas = I->deltas + (l - 1)*I->r;
    fq_zech_mpoly_struct * q = I->q + l;
    fq_zech_mpoly_struct * qt = I->qt + l;
    fq_zech_mpoly_struct * newt = I->newt + l;
    fq_zech_mpolyv_struct * delta_coeffs = I->delta_coeffs + l*I->r;

    FLINT_ASSERT(l >= 0);

    if (!fq_zech_mpoly_repack_bits_inplace(t, I->bits, ctx))
        return -1;

    if (l < 1)
    {
        for (i = 0; i < I->r; i++)
        {
            fq_zech_mpoly_divrem(I->Q, I->R, t, I->dbetas_mvar + i, ctx);
            fq_zech_mpoly_mul(I->T, I->R, I->inv_prod_dbetas_mvar + i, ctx);
            fq_zech_mpoly_divrem(I->Q, deltas + i, I->T, I->dbetas_mvar + i, ctx);
        }
        return 1;
    }

    for (i = 0; i < I->r; i++)
        delta_coeffs[i].length = 0;

    for (k = 0; k <= degs[l]; k++)
    {
        fq_zech_mpoly_divrem(q, newt, t, I->xalpha + l, ctx);
        fq_zech_mpoly_swap(t, q, ctx);
        for (j = 0; j < k; j++)
        for (i = 0; i < I->r; i++)
        {
            if (j >= delta_coeffs[i].length)
                continue;
            if (k - j >= I->prod_mbetas_coeffs[l*I->r + i].length)
                continue;

            fq_zech_mpoly_mul(qt, delta_coeffs[i].coeffs + j,
                        I->prod_mbetas_coeffs[l*I->r + i].coeffs + k - j, ctx);
            fq_zech_mpoly_sub(q, newt, qt, ctx);
            fq_zech_mpoly_swap(newt, q, ctx);
        }

        success = fq_zech_mpoly_pfrac(l - 1, newt, degs, I, ctx);
        if (success < 1)
            return success;

        for (i = 0; i < I->r; i++)
        {
            if (fq_zech_mpoly_is_zero(newdeltas + i, ctx))
                continue;

            if (k + I->prod_mbetas_coeffs[l*I->r + i].length - 1 > degs[l])
                return 0;

            fq_zech_mpolyv_set_coeff(delta_coeffs + i, k, newdeltas + i, ctx);
        }
    }

    for (i = 0; i < I->r; i++)
        fq_zech_mpoly_from_mpolyv(deltas + i, delta_coeffs + i, I->xalpha + l, ctx);

    return 1;
}
