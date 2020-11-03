/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


static void _to_polyq(
    fmpq_poly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong mask;
    slong shift, off, N;
    slong Blen = B->length;
    fmpz * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    fmpq_poly_zero(A);

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, B->bits, ctx->minfo);

    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    for (i = 0; i < Blen; i++)
        fmpq_poly_set_coeff_fmpz(A, (Bexp[N*i + off] >> shift) & mask, Bcoeff + i);
}


static int _from_polyq(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpq_poly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong var = 0;
    slong N;
    slong k;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * strideexp;
    TMP_INIT;

    if (B->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return 1;
    }

    if (!fmpz_is_one(fmpq_poly_denref(B)))
    {
        return 0;
    }

    TMP_START;

    FLINT_ASSERT(Abits <= FLINT_BITS);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    strideexp = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(strideexp, var, Abits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = B->length - 1; k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
        if (!fmpz_is_zero(B->coeffs + k))
        {
            fmpz_swap(Acoeff + Alen, B->coeffs + k);
            mpoly_monomial_mul_ui(Aexp + N*Alen, strideexp, N, k);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;

    return 1;
}


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
    fmpq_poly_t G, S, pq;

    FLINT_ASSERT(bits <= FLINT_BITS);

    I->bits = bits;
    I->r = r;
    I->w = w;

    fmpq_poly_init(I->dtq);
    fmpq_poly_init(I->S);
    fmpq_poly_init(I->R);

    I->dbetas = FLINT_ARRAY_ALLOC(r, fmpq_poly_struct);
    I->inv_prod_dbetas = FLINT_ARRAY_ALLOC(r, fmpq_poly_struct);
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

    fmpq_poly_init(G);
    fmpq_poly_init(S);
    fmpq_poly_init(pq);

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
    for (j = 0; j < r; j++)
    {
        fmpq_poly_init(I->dbetas + j);
        _to_polyq(I->dbetas + j, I->mbetas + 0*r + j, ctx);
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

    for (j = 0; j < r; j++)
        fmpq_poly_init(I->inv_prod_dbetas + j);

    for (j = 0; success && j < r; j++)
    {
        if (fmpq_poly_degree(I->dbetas + j) !=
                 fmpz_mpoly_degree_si(betas + j, 0, ctx))
        {
            success = 0;
        }
    }

    for (j = 0; success && j < r; j++)
    {
        fmpq_poly_one(pq);
        for (k = 0; k < r; k++)
        {
            if (k == j)
                continue;
            fmpq_poly_mul(pq, pq, I->dbetas + k);
        }
        fmpq_poly_xgcd(G, S, I->inv_prod_dbetas + j, I->dbetas + j, pq);
        if (!fmpq_poly_is_one(G))
            success = 0;
    }

    fmpq_poly_clear(G);
    fmpq_poly_clear(S);
    fmpq_poly_clear(pq);

    FLINT_ASSERT(success == 1);

    return success;
}


void fmpz_mpoly_pfrac_clear(
    fmpz_mpoly_pfrac_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    fmpq_poly_clear(I->dtq);
    fmpq_poly_clear(I->S);
    fmpq_poly_clear(I->R);

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
        fmpq_poly_clear(I->inv_prod_dbetas + j);
        fmpq_poly_clear(I->dbetas + j);
        for (i = 0; i <= I->w; i++)
        {
            fmpz_mpolyv_clear(I->prod_mbetas_coeffs + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->prod_mbetas + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->mbetas + i*I->r + j, ctx);
            fmpz_mpoly_clear(I->deltas + i*I->r + j, ctx);
        }
    }

    flint_free(I->inv_prod_dbetas);
    flint_free(I->dbetas);
    flint_free(I->prod_mbetas);
    flint_free(I->prod_mbetas_coeffs);
    flint_free(I->mbetas);
    flint_free(I->deltas);
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
        _to_polyq(I->dtq, t, ctx);

        success = 1;
        for (i = 0; i < I->r; i++)
        {
            fmpq_poly_mul(I->S, I->dtq, I->inv_prod_dbetas + i);
            fmpq_poly_rem(I->R, I->S, I->dbetas + i);
            if (!_from_polyq(deltas + i, I->bits, I->R, ctx))
                return 0;
        }
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

