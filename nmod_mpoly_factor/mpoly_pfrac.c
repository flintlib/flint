/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int nmod_mpoly_pfrac_init(
    nmod_mpoly_pfrac_t I,
    flint_bitcnt_t bits,
    slong r,
    slong w,
    const nmod_mpoly_struct * betas,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i, j, k;
    n_poly_t p;
    n_poly_t G, S, pq;

    FLINT_ASSERT(bits <= FLINT_BITS);

    I->bits = bits;
    I->r = r;
    I->w = w;

    I->dbetas = FLINT_ARRAY_ALLOC(r, n_poly_struct);
    I->inv_prod_dbetas = FLINT_ARRAY_ALLOC(r, n_poly_struct);
    I->prod_mbetas = FLINT_ARRAY_ALLOC((w + 1)*r, nmod_mpoly_struct);
    I->prod_mbetas_coeffs = FLINT_ARRAY_ALLOC((w + 1)*r, nmod_mpolyv_struct);
    I->mbetas = FLINT_ARRAY_ALLOC((w + 1)*r, nmod_mpoly_struct);
    I->deltas = FLINT_ARRAY_ALLOC((w + 1)*r, nmod_mpoly_struct);
    I->xalpha = FLINT_ARRAY_ALLOC(w + 1, nmod_mpoly_struct);
    I->q = FLINT_ARRAY_ALLOC(w + 1, nmod_mpoly_struct);
    I->G = FLINT_ARRAY_ALLOC(w + 1, nmod_mpoly_geobucket_struct);
    I->qt = FLINT_ARRAY_ALLOC(w + 1, nmod_mpoly_struct);
    I->newt = FLINT_ARRAY_ALLOC(w + 1, nmod_mpoly_struct);
    I->delta_coeffs = FLINT_ARRAY_ALLOC((w + 1)*r, nmod_mpolyv_struct);

    for (i = 0; i <= w; i++)
    {
        nmod_mpoly_init(I->xalpha + i, ctx);
        nmod_mpoly_init(I->q + i, ctx);
        nmod_mpoly_geobucket_init(I->G + i, ctx);
        nmod_mpoly_init(I->qt + i, ctx);
        nmod_mpoly_init(I->newt + i, ctx);
        for (j = 0; j < r; j++)
        {
            nmod_mpoly_init(I->deltas + i*r + j, ctx);
            nmod_mpolyv_init(I->delta_coeffs + i*r + j, ctx);
        }

        if (i < 1)
            continue;

        nmod_mpoly_gen(I->xalpha + i, i, ctx);
        nmod_mpoly_sub_ui(I->xalpha + i, I->xalpha + i, alpha[i - 1], ctx);
        nmod_mpoly_repack_bits_inplace(I->xalpha + i, I->bits, ctx);
    }

    n_poly_init(p);
    n_poly_init(G);
    n_poly_init(S);
    n_poly_init(pq);

    /* set betas */
    i = w;
    for (j = 0; j < r; j++)
    {
        nmod_mpoly_init(I->mbetas + i*r + j, ctx);
        nmod_mpoly_set(I->mbetas + i*r + j, betas + j, ctx);
    }
    for (i--; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            nmod_mpoly_init(I->mbetas + i*r + j, ctx);
            nmod_mpoly_evaluate_one_ui(I->mbetas + i*r + j,
                              I->mbetas + (i + 1)*r + j, i + 1, alpha[i], ctx);
        }
    }

    for (j = 0; j < r; j++)
    {
        n_poly_init(I->dbetas + j);
        if (!nmod_mpoly_get_n_poly(I->dbetas + j,
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
            nmod_mpoly_init(I->prod_mbetas + i*r + j, ctx);
            nmod_mpoly_one(I->prod_mbetas + i*r + j, ctx);
            for (k = 0; k < r; k++)
            {
                if (k == j)
                    continue;
                nmod_mpoly_mul(I->prod_mbetas + i*r + j,
                           I->prod_mbetas + i*r + j, I->mbetas + i*r + k, ctx);
            }
            nmod_mpolyv_init(I->prod_mbetas_coeffs + i*r + j, ctx);
            if (i > 0)
            {
                nmod_mpoly_to_mpolyv(I->prod_mbetas_coeffs + i*r + j,
                                 I->prod_mbetas + i*r + j, I->xalpha + i, ctx);
            }
        }
    }

    for (j = 0; j < r; j++)
        n_poly_init(I->inv_prod_dbetas + j);

    for (j = 0; success && j < r; j++)
    {
        if (n_poly_degree(I->dbetas + j) !=
            nmod_mpoly_degree_si(betas + j, 0, ctx))
        {
            success = 0;
        }
    }

    for (j = 0; success && j < r; j++)
    {
        n_poly_one(pq);
        for (k = 0; k < r; k++)
        {
            if (k == j)
                continue;
            n_poly_mod_mul(pq, pq, I->dbetas + k, ctx->mod);
        }
        n_poly_mod_xgcd(G, S, I->inv_prod_dbetas + j, I->dbetas + j, pq, ctx->mod);
        if (!n_poly_is_one(G))
        {
            success = 0;
        }
    }

    n_poly_clear(p);
    n_poly_clear(G);
    n_poly_clear(S);
    n_poly_clear(pq);

    I->dbetas_mvar = FLINT_ARRAY_ALLOC(r, nmod_mpoly_struct);
    I->inv_prod_dbetas_mvar = FLINT_ARRAY_ALLOC(r, nmod_mpoly_struct);
    for (j = 0; j < r; j++)
    {
        nmod_mpoly_init(I->dbetas_mvar + j, ctx);
        nmod_mpoly_init(I->inv_prod_dbetas_mvar + j, ctx);

        _nmod_mpoly_set_nmod_poly(I->dbetas_mvar + j, I->bits,
                             I->dbetas[j].coeffs, I->dbetas[j].length, 0, ctx);

        _nmod_mpoly_set_nmod_poly(I->inv_prod_dbetas_mvar + j, I->bits,
           I->inv_prod_dbetas[j].coeffs, I->inv_prod_dbetas[j].length, 0, ctx);
    }

    nmod_mpoly_init(I->T, ctx);
    nmod_mpoly_init(I->Q, ctx);
    nmod_mpoly_init(I->R, ctx);

    FLINT_ASSERT(success == 1);

    return success;
}

void nmod_mpoly_pfrac_clear(
    nmod_mpoly_pfrac_t I,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i <= I->w; i++)
    {
        nmod_mpoly_clear(I->xalpha + i, ctx);
        nmod_mpoly_clear(I->q + i, ctx);
        nmod_mpoly_geobucket_clear(I->G + i, ctx);
        nmod_mpoly_clear(I->qt + i, ctx);
        nmod_mpoly_clear(I->newt + i, ctx);
        for (j = 0; j < I->r; j++)
            nmod_mpolyv_clear(I->delta_coeffs + i*I->r + j, ctx);
    }

    flint_free(I->xalpha);
    flint_free(I->q);
    flint_free(I->G);
    flint_free(I->qt);
    flint_free(I->newt);
    flint_free(I->delta_coeffs);

    for (j = 0; j < I->r; j++)
    {
        n_poly_clear(I->inv_prod_dbetas + j);
        n_poly_clear(I->dbetas + j);
        for (i = 0; i <= I->w; i++)
        {
            nmod_mpolyv_clear(I->prod_mbetas_coeffs + i*I->r + j, ctx);
            nmod_mpoly_clear(I->prod_mbetas + i*I->r + j, ctx);
            nmod_mpoly_clear(I->mbetas + i*I->r + j, ctx);
            nmod_mpoly_clear(I->deltas + i*I->r + j, ctx);
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
        nmod_mpoly_clear(I->dbetas_mvar + j, ctx);
        nmod_mpoly_clear(I->inv_prod_dbetas_mvar + j, ctx);
    }
    flint_free(I->dbetas_mvar);
    flint_free(I->inv_prod_dbetas_mvar);

    nmod_mpoly_clear(I->T, ctx);
    nmod_mpoly_clear(I->Q, ctx);
    nmod_mpoly_clear(I->R, ctx);
}


int nmod_mpoly_pfrac(
    slong l,
    nmod_mpoly_t t,
    const slong * degs,
    nmod_mpoly_pfrac_t I,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    nmod_mpoly_struct * deltas = I->deltas + l*I->r;
    nmod_mpoly_struct * newdeltas = I->deltas + (l - 1)*I->r;
    nmod_mpoly_struct * q = I->q + l;
    nmod_mpoly_struct * qt = I->qt + l;
    nmod_mpoly_struct * newt = I->newt + l;
    nmod_mpolyv_struct * delta_coeffs = I->delta_coeffs + l*I->r;
    nmod_mpoly_geobucket_struct * G = I->G + l;

    FLINT_ASSERT(l >= 0);

    if (!nmod_mpoly_repack_bits_inplace(t, I->bits, ctx))
        return -1;

    if (l < 1)
    {
        for (i = 0; i < I->r; i++)
        {
            nmod_mpoly_divrem(I->Q, I->R, t, I->dbetas_mvar + i, ctx);
            nmod_mpoly_mul(I->T, I->R, I->inv_prod_dbetas_mvar + i, ctx);
            nmod_mpoly_divrem(I->Q, deltas + i, I->T, I->dbetas_mvar + i, ctx);
        }
        return 1;
    }

    for (i = 0; i < I->r; i++)
        delta_coeffs[i].length = 0;

    for (k = 0; k <= degs[l]; k++)
    {
        nmod_mpoly_divrem(q, newt, t, I->xalpha + l, ctx);
        nmod_mpoly_swap(t, q, ctx);
        nmod_mpoly_geobucket_set(G, newt, ctx);

        for (j = 0; j < k; j++)
        for (i = 0; i < I->r; i++)
        {
            if (j >= delta_coeffs[i].length)
                continue;
            if (k - j >= I->prod_mbetas_coeffs[l*I->r + i].length)
                continue;

            nmod_mpoly_mul(qt, delta_coeffs[i].coeffs + j,
                        I->prod_mbetas_coeffs[l*I->r + i].coeffs + k - j, ctx);
            nmod_mpoly_geobucket_sub(G, qt, ctx);
        }

        nmod_mpoly_geobucket_empty(newt, G, ctx);

        if (nmod_mpoly_is_zero(newt, ctx))
            continue;

        success = nmod_mpoly_pfrac(l - 1, newt, degs, I, ctx);
        if (success < 1)
            return success;

        for (i = 0; i < I->r; i++)
        {
            if (nmod_mpoly_is_zero(newdeltas + i, ctx))
                continue;

            if (k + I->prod_mbetas_coeffs[l*I->r + i].length - 1 > degs[l])
                return 0;

            nmod_mpolyv_set_coeff(delta_coeffs + i, k, newdeltas + i, ctx);
        }
    }

    for (i = 0; i < I->r; i++)
        nmod_mpoly_from_mpolyv(deltas + i, I->bits,
                                         delta_coeffs + i, I->xalpha + l, ctx);

    return 1;
}
