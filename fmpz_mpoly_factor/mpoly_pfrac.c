/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_pfrac_init(
    fmpz_mpoly_pfrac_t I,
    flint_bitcnt_t bits,
    slong r,
    slong w,
    const fmpz_mpoly_struct * betas,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong success = 1;
    slong i, j, k;

    FLINT_ASSERT(bits <= FLINT_BITS);

    I->bits = bits;
    I->r = r;
    I->w = w;

    I->prod_mbetas = FLINT_ARRAY_ALLOC((w + 1)*r, fmpz_mpoly_struct);
    I->prod_mbetas_coeffs = FLINT_ARRAY_ALLOC((w + 1)*r, fmpz_mpolyv_struct);
    I->mbetas = FLINT_ARRAY_ALLOC((w + 1)*r, fmpz_mpoly_struct);
    I->deltas = FLINT_ARRAY_ALLOC((w + 1)*r, fmpz_mpoly_struct);
    I->xalpha = FLINT_ARRAY_ALLOC(w + 1, fmpz_mpoly_struct);
    I->q = FLINT_ARRAY_ALLOC(w + 1, fmpz_mpoly_struct);
    I->U = FLINT_ARRAY_ALLOC(w + 1, fmpz_mpoly_univar_struct);
    I->G = FLINT_ARRAY_ALLOC(w + 1, fmpz_mpoly_geobucket_struct);
    I->qt = FLINT_ARRAY_ALLOC(w + 1, fmpz_mpoly_struct);
    I->newt = FLINT_ARRAY_ALLOC(w + 1, fmpz_mpoly_struct);
    I->delta_coeffs = FLINT_ARRAY_ALLOC((w + 1)*r, fmpz_mpolyv_struct);

    for (i = 0; i <= w; i++)
    {
        fmpz_mpoly_init(I->xalpha + i, ctx);
        fmpz_mpoly_init(I->q + i, ctx);
        fmpz_mpoly_univar_init(I->U + i, ctx);
        fmpz_mpoly_geobucket_init(I->G + i, ctx);
        fmpz_mpoly_init(I->qt + i, ctx);
        fmpz_mpoly_init(I->newt + i, ctx);
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_init(I->deltas + i*r + j, ctx);
            fmpz_mpolyv_init(I->delta_coeffs + i*r + j, ctx);
        }

        if (i < 1)
            continue;

        fmpz_mpoly_gen(I->xalpha + i, i, ctx);
        fmpz_mpoly_sub_fmpz(I->xalpha + i, I->xalpha + i, alpha + i - 1, ctx);
        fmpz_mpoly_repack_bits_inplace(I->xalpha + i, I->bits, ctx);
    }

    /* set betas */
    i = w;
    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_init(I->mbetas + i*r + j, ctx);
        fmpz_mpoly_set(I->mbetas + i*r + j, betas + j, ctx);
    }
    for (i--; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_init(I->mbetas + i*r + j, ctx);
            fmpz_mpoly_evaluate_one_fmpz(I->mbetas + i*r + j,
                             I->mbetas + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    /* set product of betas */
    for (i = w; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_init(I->prod_mbetas + i*r + j, ctx);
            fmpz_mpoly_one(I->prod_mbetas + i*r + j, ctx);
            for (k = 0; k < r; k++)
            {
                if (k == j)
                    continue;
                fmpz_mpoly_mul(I->prod_mbetas + i*r + j,
                           I->prod_mbetas + i*r + j, I->mbetas + i*r + k, ctx);
            }
            fmpz_mpolyv_init(I->prod_mbetas_coeffs + i*r + j, ctx);
            if (i > 0)
            {
                fmpz_mpoly_to_mpolyv(I->prod_mbetas_coeffs + i*r + j,
                                 I->prod_mbetas + i*r + j, I->xalpha + i, ctx);
            }
        }        
    }

    fmpz_poly_pfrac_init(I->uni_pfrac);
    fmpz_poly_init(I->uni_a);
    I->uni_c = FLINT_ARRAY_ALLOC(r, fmpz_poly_struct);
    for (j = 0; j < r; j++)
    {
        fmpz_poly_init(I->uni_c + j);
        fmpz_mpoly_get_fmpz_poly(I->uni_c + j, I->mbetas + 0*r + j, 0, ctx);

        success = success && (fmpz_poly_degree(I->uni_c + j) ==
                                      fmpz_mpoly_degree_si(betas + j, 0, ctx));
    }

    success = success && fmpz_poly_pfrac_precompute(I->uni_pfrac, I->uni_c, r);

    if (!success)
        flint_throw(FLINT_ERROR, "fmpz_mpoly_pfrac_init: internal error");

    return success;
}


void fmpz_mpoly_pfrac_clear(
    fmpz_mpoly_pfrac_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i <= I->w; i++)
    {
        fmpz_mpoly_clear(I->xalpha + i, ctx);
        fmpz_mpoly_clear(I->q + i, ctx);
        fmpz_mpoly_univar_clear(I->U + i, ctx);
        fmpz_mpoly_geobucket_clear(I->G + i, ctx);
        fmpz_mpoly_clear(I->qt + i, ctx);
        fmpz_mpoly_clear(I->newt + i, ctx);
        for (j = 0; j < I->r; j++)
            fmpz_mpolyv_clear(I->delta_coeffs + i*I->r + j, ctx);
    }
    flint_free(I->xalpha);
    flint_free(I->q);
    flint_free(I->U);
    flint_free(I->G);
    flint_free(I->qt);
    flint_free(I->newt);
    flint_free(I->delta_coeffs);

    for (j = 0; j < I->r; j++)
    {
        for (i = 0; i <= I->w; i++)
        {
            fmpz_mpolyv_clear(I->prod_mbetas_coeffs + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->prod_mbetas + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->mbetas + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->deltas + i*I->r + j, ctx);
        }
    }

    flint_free(I->prod_mbetas);
    flint_free(I->prod_mbetas_coeffs);
    flint_free(I->mbetas);
    flint_free(I->deltas);

    fmpz_poly_pfrac_clear(I->uni_pfrac);
    fmpz_poly_clear(I->uni_a);
    for (j = 0; j < I->r; j++)
        fmpz_poly_clear(I->uni_c + j);
    flint_free(I->uni_c);
}


int fmpz_mpoly_pfrac(
    slong l,
    fmpz_mpoly_t t,
    const slong * degs,
    fmpz_mpoly_pfrac_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    int success, use_U;
    slong i, j, k, Ui;
    fmpz_mpoly_struct * deltas = I->deltas + l*I->r;
    fmpz_mpoly_struct * newdeltas = I->deltas + (l - 1)*I->r;
    fmpz_mpoly_struct * q = I->q + l;
    fmpz_mpoly_struct * qt = I->qt + l;
    fmpz_mpoly_struct * newt = I->newt + l;
    fmpz_mpolyv_struct * delta_coeffs = I->delta_coeffs + l*I->r;
    fmpz_mpoly_geobucket_struct * G = I->G + l;
    fmpz_mpoly_univar_struct * U = I->U + l;

    FLINT_ASSERT(l >= 0);

    if (!fmpz_mpoly_repack_bits_inplace(t, I->bits, ctx))
        return -1;

    if (l < 1)
    {
        fmpz_mpoly_get_fmpz_poly(I->uni_a, t, 0, ctx);

        if (!fmpz_poly_pfrac_precomp(I->uni_c, I->uni_a, I->uni_pfrac))
            return 0;

        for (i = 0; i < I->r; i++)
            _fmpz_mpoly_set_fmpz_poly(deltas + i, I->bits,
                               I->uni_c[i].coeffs, I->uni_c[i].length, 0, ctx);
        return 1;
    }

    for (i = 0; i < I->r; i++)
        delta_coeffs[i].length = 0;

    use_U = I->xalpha[l].length == 1;
    if (use_U)
        fmpz_mpoly_to_univar(U, t, l, ctx);
    Ui = U->length - 1;

    for (k = 0; k <= degs[l]; k++)
    {
        if (use_U)
        {
            if (Ui >= 0 && fmpz_equal_si(U->exps + Ui, k))
            {
                fmpz_mpoly_geobucket_set(G, U->coeffs + Ui, ctx);
                Ui--;
            }
            else
            {
                G->length = 0;
            }
        }
        else
        {
            fmpz_mpoly_divrem(q, newt, t, I->xalpha + l, ctx);
            fmpz_mpoly_swap(t, q, ctx);
            fmpz_mpoly_geobucket_set(G, newt, ctx);
        }

        for (j = 0; j < k; j++)
        for (i = 0; i < I->r; i++)
        {
            if (j >= delta_coeffs[i].length)
                continue;
            if (k - j >= I->prod_mbetas_coeffs[l*I->r + i].length)
                continue;

            fmpz_mpoly_mul(qt, delta_coeffs[i].coeffs + j,
                        I->prod_mbetas_coeffs[l*I->r + i].coeffs + k - j, ctx);
            fmpz_mpoly_geobucket_sub(G, qt, ctx);
        }

        fmpz_mpoly_geobucket_empty(newt, G, ctx);

        if (fmpz_mpoly_is_zero(newt, ctx))
            continue;

        success = fmpz_mpoly_pfrac(l - 1, newt, degs, I, ctx);
        if (success < 1)
            return success;

        for (i = 0; i < I->r; i++)
        {
            if (fmpz_mpoly_is_zero(newdeltas + i, ctx))
                continue;

            if (k + I->prod_mbetas_coeffs[l*I->r + i].length - 1 > degs[l])
                return 0;

            fmpz_mpolyv_set_coeff(delta_coeffs + i, k, newdeltas + i, ctx);
        }
    }

    for (i = 0; i < I->r; i++)
        fmpz_mpoly_from_mpolyv(deltas + i, I->bits,
                                         delta_coeffs + i, I->xalpha + l, ctx);

    return 1;
}

