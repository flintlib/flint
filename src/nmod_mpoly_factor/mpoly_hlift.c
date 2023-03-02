/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

static int _hlift_quartic2(
    slong m,
    nmod_mpoly_struct * f,
    slong r,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    nmod_mpoly_t Aq, t, t2, t3, xalpha;
    nmod_mpoly_geobucket_t G;
    nmod_mpoly_struct betas[2], * deltas;
    nmod_mpoly_pfrac_t I;
    nmod_mpolyv_struct B[2];
    slong tdeg;
    flint_bitcnt_t bits = A->bits;

    FLINT_ASSERT(r == 2);
    r = 2;

    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(t2, ctx);
    nmod_mpoly_init(t3, ctx);
    nmod_mpoly_init(xalpha, ctx);
    nmod_mpoly_init(Aq, ctx);
    nmod_mpoly_geobucket_init(G, ctx);

    nmod_mpoly_gen(xalpha, m, ctx);
    nmod_mpoly_sub_ui(xalpha, xalpha, alpha[m - 1], ctx);
    nmod_mpoly_repack_bits_inplace(xalpha, bits, ctx);

    for (i = 0; i < r; i++)
    {
        nmod_mpolyv_init(B + i, ctx);
        nmod_mpoly_repack_bits_inplace(f + i, bits, ctx);
        nmod_mpoly_to_mpolyv(B + i, f + i, xalpha, ctx);
        nmod_mpolyv_fit_length(B + i, degs[m] + 1, ctx);
        for (j = B[i].length; j <= degs[m]; j++)
            nmod_mpoly_zero(B[i].coeffs + j, ctx);
    }

    for (i = 0; i < r; i++)
        betas[i] = B[i].coeffs[0];

    success = nmod_mpoly_pfrac_init(I, bits, r, m - 1, betas, alpha, ctx);
    FLINT_ASSERT(success == 1);
    deltas = I->deltas + (m - 1)*I->r;

    nmod_mpoly_divrem(Aq, t, A, xalpha, ctx);

#if FLINT_WANT_ASSERT
    nmod_mpoly_one(t2, ctx);
    for (i = 0; i < r; i++)
        nmod_mpoly_mul(t2, t2, betas + i, ctx);
    FLINT_ASSERT(nmod_mpoly_equal(t, t2, ctx));
#endif

    for (j = 1; j <= degs[m]; j++)
    {
        nmod_mpoly_divrem(t2, t, Aq, xalpha, ctx);
        nmod_mpoly_swap(Aq, t2, ctx);
        nmod_mpoly_geobucket_set(G, t, ctx);

        for (i = 0; i <= j; i++)
        {
            nmod_mpoly_mul(t, B[0].coeffs + i, B[1].coeffs + j - i, ctx);
            nmod_mpoly_geobucket_sub(G, t, ctx);
        }
        nmod_mpoly_geobucket_empty(t, G, ctx);

        if (nmod_mpoly_is_zero(t, ctx))
            continue;

        success = nmod_mpoly_pfrac(m - 1, t, degs, I, ctx);
        if (success < 1)
        {
            success = 0;
            goto cleanup;
        }

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            nmod_mpoly_add(t3, B[i].coeffs + j, deltas + i, ctx);
            nmod_mpoly_swap(B[i].coeffs + j, t3, ctx);
            if (!nmod_mpoly_is_zero(B[i].coeffs + j, ctx))
                B[i].length = FLINT_MAX(B[i].length, j + 1);
            FLINT_ASSERT(B[i].length > 0);
            tdeg += B[i].length - 1;
        }

        if (tdeg > degs[m])
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    nmod_mpoly_pfrac_clear(I, ctx);

    for (i = 0; i < r; i++)
    {
        if (success)
            nmod_mpoly_from_mpolyv(f + i, bits, B + i, xalpha, ctx);
        nmod_mpolyv_clear(B + i, ctx);
    }

    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_clear(t2, ctx);
    nmod_mpoly_clear(t3, ctx);
    nmod_mpoly_clear(xalpha, ctx);
    nmod_mpoly_clear(Aq, ctx);
    nmod_mpoly_geobucket_clear(G, ctx);

    return success;
}


static int _hlift_quartic(
    slong m,
    nmod_mpoly_struct * f,
    slong r,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    nmod_mpoly_t Aq, t, t1, t2, t3, xalpha;
    nmod_mpoly_geobucket_t G;
    nmod_mpoly_struct * betas, * deltas;
    nmod_mpoly_pfrac_t I;
    nmod_mpolyv_struct * B, * U;
    slong tdeg;
    flint_bitcnt_t bits = A->bits;

    FLINT_ASSERT(r > 2);

    B = FLINT_ARRAY_ALLOC(2*r, nmod_mpolyv_struct);
    U = B + r;

    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);
    nmod_mpoly_init(t3, ctx);
    nmod_mpoly_init(xalpha, ctx);
    nmod_mpoly_init(Aq, ctx);
    nmod_mpoly_geobucket_init(G, ctx);

    nmod_mpoly_gen(xalpha, m, ctx);
    nmod_mpoly_sub_ui(xalpha, xalpha, alpha[m - 1], ctx);
    nmod_mpoly_repack_bits_inplace(xalpha, bits, ctx);

    for (k = 0; k < r; k++)
    {
        nmod_mpolyv_init(U + k, ctx);
        nmod_mpolyv_fit_length(U + k, degs[m] + 1, ctx);
        for (j = 0; j <= degs[m]; j++)
            nmod_mpoly_zero(U[k].coeffs + j, ctx);

        nmod_mpolyv_init(B + k, ctx);
        nmod_mpoly_repack_bits_inplace(f + k, bits, ctx);
        nmod_mpoly_to_mpolyv(B + k, f + k, xalpha, ctx);
        nmod_mpolyv_fit_length(B + k, degs[m] + 1, ctx);
        for (j = B[k].length; j <= degs[m]; j++)
            nmod_mpoly_zero(B[k].coeffs + j, ctx);
    }

    betas = FLINT_ARRAY_ALLOC(r, nmod_mpoly_struct);
    for (i = 0; i < r; i++)
        betas[i] = B[i].coeffs[0];

    success = nmod_mpoly_pfrac_init(I, bits, r, m - 1, betas, alpha, ctx);
    FLINT_ASSERT(success == 1);
    deltas = I->deltas + (m - 1)*I->r;

    k = r - 2;
    nmod_mpoly_mul(U[k].coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, ctx);
    for (k--; k >= 1; k--)
        nmod_mpoly_mul(U[k].coeffs + 0, B[k].coeffs + 0, U[k + 1].coeffs + 0, ctx);

    nmod_mpoly_divrem(t2, t, A, xalpha, ctx);
    nmod_mpoly_swap(Aq, t2, ctx);

#if FLINT_WANT_ASSERT
    nmod_mpoly_one(t2, ctx);
    for (i = 0; i < r; i++)
        nmod_mpoly_mul(t2, t2, betas + i, ctx);
    FLINT_ASSERT(nmod_mpoly_equal(t, t2, ctx));
#endif

    for (j = 1; j <= degs[m]; j++)
    {
        k = r - 2;

        G->length = 0;
        for (i = 0; i <= j; i++)
        {
            nmod_mpoly_mul(t1, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
            nmod_mpoly_geobucket_add(G, t1, ctx);
        }
        nmod_mpoly_geobucket_empty(U[k].coeffs + j, G, ctx);

        for (k--; k >= 1; k--)
        {
            G->length = 0;
            for (i = 0; i <= j; i++)
            {
                nmod_mpoly_mul(t1, B[k].coeffs + i, U[k + 1].coeffs + j - i, ctx);
                nmod_mpoly_geobucket_add(G, t1, ctx);
            }
            nmod_mpoly_geobucket_empty(U[k].coeffs + j, G, ctx);
        }

        nmod_mpoly_divrem(t2, t, Aq, xalpha, ctx);
        nmod_mpoly_swap(Aq, t2, ctx);
        nmod_mpoly_geobucket_set(G, t, ctx);

        for (i = 0; i <= j; i++)
        {
            nmod_mpoly_mul(t, B[0].coeffs + i, U[1].coeffs + j - i, ctx);
            nmod_mpoly_geobucket_sub(G, t, ctx);
        }
        nmod_mpoly_geobucket_empty(t, G, ctx);

        if (nmod_mpoly_is_zero(t, ctx))
            continue;

        success = nmod_mpoly_pfrac(m - 1, t, degs, I, ctx);
        if (success < 1)
        {
            success = 0;
            goto cleanup;
        }

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            nmod_mpoly_add(t3, B[i].coeffs + j, deltas + i, ctx);
            nmod_mpoly_swap(B[i].coeffs + j, t3, ctx);
            if (!nmod_mpoly_is_zero(B[i].coeffs + j, ctx))
                B[i].length = FLINT_MAX(B[i].length, j + 1);
            FLINT_ASSERT(B[i].length > 0);
            tdeg += B[i].length - 1;
        }

        if (tdeg > degs[m])
        {
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        nmod_mpoly_mul(t, B[k].coeffs + 0, deltas + k + 1, ctx);
        nmod_mpoly_mul(t1, deltas + k, B[k + 1].coeffs + 0, ctx);
        nmod_mpoly_add(t, t, t1, ctx);
        nmod_mpoly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
        for (k--; k >= 1; k--)
        {
            nmod_mpoly_mul(t1, B[k].coeffs + 0, t, ctx);
            nmod_mpoly_swap(t, t1, ctx);
            nmod_mpoly_mul(t1, deltas + k, U[k + 1].coeffs + 0, ctx);
            nmod_mpoly_add(t, t, t1, ctx);
            nmod_mpoly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
        }
    }

    success = 1;

cleanup:

    nmod_mpoly_pfrac_clear(I, ctx);

    flint_free(betas);

    for (i = 0; i < r; i++)
    {
        if (success)
            nmod_mpoly_from_mpolyv(f + i, bits, B + i, xalpha, ctx);
        nmod_mpolyv_clear(B + i, ctx);
        nmod_mpolyv_clear(U + i, ctx);
    }

    flint_free(B);

    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_clear(t1, ctx);
    nmod_mpoly_clear(t2, ctx);
    nmod_mpoly_clear(t3, ctx);
    nmod_mpoly_clear(xalpha, ctx);
    nmod_mpoly_clear(Aq, ctx);
    nmod_mpoly_geobucket_clear(G, ctx);

    return success;
}


static int _hlift_quintic(
    slong m,
    nmod_mpoly_struct * f,
    slong r,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    nmod_mpoly_t e, t, pow, xalpha, q;
    nmod_mpoly_struct * betas, * deltas;
    nmod_mpoly_pfrac_t I;
    flint_bitcnt_t bits = A->bits;

    FLINT_ASSERT(r > 1);

    nmod_mpoly_init(e, ctx);
    nmod_mpoly_init(t, ctx);
    nmod_mpoly_init(pow, ctx);
    nmod_mpoly_init(xalpha, ctx);
    nmod_mpoly_init(q, ctx);

    betas = FLINT_ARRAY_ALLOC(r, nmod_mpoly_struct);
    for (i = 0; i < r; i++)
    {
        nmod_mpoly_init(betas + i, ctx);
        nmod_mpoly_repack_bits_inplace(f + i, bits, ctx);
        nmod_mpoly_evaluate_one_ui(betas + i, f + i, m, alpha[m - 1], ctx);
    }

    nmod_mpoly_mul(t, f + 0, f + 1, ctx);
    for (i = 2; i < r; i++)
        nmod_mpoly_mul(t, t, f + i, ctx);
    nmod_mpoly_sub(e, A, t, ctx);

    nmod_mpoly_one(pow, ctx);
    nmod_mpoly_repack_bits_inplace(pow, bits, ctx);

    nmod_mpoly_gen(xalpha, m, ctx);
    nmod_mpoly_sub_ui(xalpha, xalpha, alpha[m - 1], ctx);
    nmod_mpoly_repack_bits_inplace(xalpha, bits, ctx);

    nmod_mpoly_pfrac_init(I, bits, r, m - 1, betas, alpha, ctx);
    deltas = I->deltas + (m - 1)*I->r;

    for (j = 1; j <= degs[m]; j++)
    {
        if (nmod_mpoly_is_zero(e, ctx))
        {
            success = 1;
            goto cleanup;
        }

        nmod_mpoly_mul(pow, pow, xalpha, ctx);
        success = nmod_mpoly_divides(q, e, pow, ctx);
        FLINT_ASSERT(success);
        nmod_mpoly_evaluate_one_ui(t, q, m, alpha[m - 1], ctx);

        success = nmod_mpoly_pfrac(m - 1, t, degs, I, ctx);
        if (success < 1)
        {
            success = 0;
            goto cleanup;
        }

        for (i = 0; i < r; i++)
        {
            nmod_mpoly_mul(t, deltas + i, pow, ctx);
            nmod_mpoly_add(f + i, f + i, t, ctx);
        }

        nmod_mpoly_mul(t, f + 0, f + 1, ctx);
        for (i = 2; i < r; i++)
            nmod_mpoly_mul(t, t, f + i, ctx);
        nmod_mpoly_sub(e, A, t, ctx);
    }

    success = nmod_mpoly_is_zero(e, ctx);

cleanup:

    nmod_mpoly_pfrac_clear(I, ctx);

    nmod_mpoly_clear(e, ctx);
    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_clear(pow, ctx);
    nmod_mpoly_clear(xalpha, ctx);
    nmod_mpoly_clear(q, ctx);

    for (i = 0; i < r; i++)
    {
        if (success)
            nmod_mpoly_repack_bits_inplace(f + i, bits, ctx);

        nmod_mpoly_clear(betas + i, ctx);
    }

    flint_free(betas);

    return success;
}

/* should have A = prod_i f[i] mod (gen(m) - alpha[m-1]) */
int nmod_mpoly_hlift(
    slong m,
    nmod_mpoly_struct * f, /* length r */
    slong r,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(r >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    if (r == 2)
        return _hlift_quartic2(m, f, r, alpha, A, degs, ctx);
    else if (r < 20)
        return _hlift_quartic(m, f, r, alpha, A, degs, ctx);
    else
        return _hlift_quintic(m, f, r, alpha, A, degs, ctx);
}
