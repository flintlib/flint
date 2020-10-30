/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"

/*
    input A, B0, B1 with A(y,x) = B0(y,x) * B1(y,x) mod (y-alpha)
    return
       -1: B0(alpha,x) & B1(alpha,x) are not pairwise prime, or
           A(alpha,x) has wrong degree w.r.t x
        0: lift of B0 and B1 to true factors is impossible
        1: successfully lifted B0 and B1 to true factors without changing lc_x
*/

int n_fq_bpoly_hlift2_cubic(
    n_fq_bpoly_t A, /* clobbered (shifted by alpha) */
    n_fq_bpoly_t B0,
    n_fq_bpoly_t B1,
    const fq_nmod_t alpha_,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx,
    nmod_eval_interp_t E,
    n_poly_bpoly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    int success;
    slong len = nmod_eval_interp_eval_length(E);
    slong i, j;
    n_fq_poly_struct * c, * s, * t, * u, * v, * g, * ce;
    n_fq_bpoly_struct * B0e, * B1e;
    mp_limb_t * alpha;

    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(B0, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(B1, ctx));

    if (A->length < 1 || B0->length < 1 || B1->length < 1)
        return -1;

    n_poly_stack_fit_request(St->poly_stack, 7);
    c  = n_poly_stack_take_top(St->poly_stack);
    s  = n_poly_stack_take_top(St->poly_stack);
    t  = n_poly_stack_take_top(St->poly_stack);
    u  = n_poly_stack_take_top(St->poly_stack);
    v  = n_poly_stack_take_top(St->poly_stack);
    g  = n_poly_stack_take_top(St->poly_stack);
    ce = n_poly_stack_take_top(St->poly_stack);

    n_bpoly_stack_fit_request(St->bpoly_stack, 2);
    B0e  = n_bpoly_stack_take_top(St->bpoly_stack);
    B1e  = n_bpoly_stack_take_top(St->bpoly_stack);

    alpha = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    n_fq_set_fq_nmod(alpha, alpha_, ctx);

    n_fq_bpoly_taylor_shift_gen0_n_fq(A, alpha, ctx);
    n_fq_bpoly_taylor_shift_gen0_n_fq(B0, alpha, ctx);
    n_fq_bpoly_taylor_shift_gen0_n_fq(B1, alpha, ctx);

#if FLINT_WANT_ASSERT
    {
        n_fq_poly_t T;
        n_fq_poly_init(T);
        n_fq_poly_mul(T, B0->coeffs + 0, B1->coeffs + 0, ctx);
        FLINT_ASSERT(n_fq_poly_equal(A->coeffs + 0, T, ctx));
        n_fq_poly_clear(T);
    }
#endif

    if (n_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(n_bpoly_degree1(A) == n_poly_degree(A->coeffs + 0));

    n_fq_poly_xgcd(g, s, t, B1->coeffs + 0, B0->coeffs + 0, ctx);
    if (!n_fq_poly_is_one(g, ctx))
    {
        success = -1;
        goto cleanup;
    }

    n_fq_bpoly_fit_length(B0, A->length);
    n_fq_bpoly_fit_length(B0e, A->length);
    for (i = 0; i < B0->length; i++)
        nmod_eval_interp_from_coeffs_n_fq_poly(B0e->coeffs + i, B0->coeffs + i, E, ctx);
    for (i = B0->length; i < A->length; i++)
    {
        n_fq_poly_zero(B0->coeffs + i);
        n_fq_evals_zero(B0e->coeffs + i);
    }

    n_fq_bpoly_fit_length(B1, A->length);
    n_fq_bpoly_fit_length(B1e, A->length);
    for (i = 0; i < B1->length; i++)
        nmod_eval_interp_from_coeffs_n_fq_poly(B1e->coeffs + i, B1->coeffs + i, E, ctx);
    for (i = B1->length; i < A->length; i++)
    {
        n_fq_poly_zero(B1->coeffs + i);
        n_fq_evals_zero(B1e->coeffs + i);
    }

    for (j = 1; j < A->length; j++)
    {
        n_fq_evals_zero(ce);
        for (i = 0; i <= j; i++)
        {
            if (i < B0->length && j - i < B1->length)
            {
                n_fq_evals_addmul(ce, B0e->coeffs + i,
                                      B1e->coeffs + j - i, len, ctx);
            }
        }

        nmod_eval_interp_to_coeffs_n_fq_poly(c, ce, E, ctx);
        n_fq_poly_sub(c, A->coeffs + j, c, ctx);

#if FLINT_WANT_ASSERT
        {
            n_fq_poly_t c_check;
            n_fq_poly_init(c_check);
            n_fq_poly_set(c_check, A->coeffs + j, ctx);
            for (i = 0; i <= j; i++)
            {
                if (i < B0->length && j - i < B1->length)
                {

                    n_fq_poly_mul(t, B0->coeffs + i, B1->coeffs + j - i, ctx);
                    n_fq_poly_sub(c_check, c_check, t, ctx);
                }
            }

            assert(n_fq_poly_equal(c, c_check, ctx));

            n_fq_poly_clear(c_check);
        }
#endif

        if (n_fq_poly_is_zero(c))
            continue;

        n_fq_poly_mul_(t, s, c, ctx, St->poly_stack);
        n_fq_poly_divrem_(g, u, t, B0->coeffs + 0, ctx, St->poly_stack);
        n_fq_poly_mul_(t, u, B1->coeffs + 0, ctx, St->poly_stack);
        n_fq_poly_sub(ce, c, t, ctx);
        n_fq_poly_divrem_(v, g, ce, B0->coeffs + 0, ctx, St->poly_stack);

        if (!n_fq_poly_is_zero(u))
        {
            n_fq_poly_add(B0->coeffs + j, B0->coeffs + j, u, ctx);
            nmod_eval_interp_from_coeffs_n_fq_poly(B0e->coeffs + j, B0->coeffs + j, E, ctx);
        }

        if (!n_fq_poly_is_zero(v))
        {
            n_fq_poly_add(B1->coeffs + j, B1->coeffs + j, v, ctx);
            nmod_eval_interp_from_coeffs_n_fq_poly(B1e->coeffs + j, B1->coeffs + j, E, ctx);
        }

        if (!n_fq_poly_is_zero(B0->coeffs + j))
            B0->length = FLINT_MAX(B0->length, j + 1);

        if (!n_fq_poly_is_zero(B1->coeffs + j))
            B1->length = FLINT_MAX(B1->length, j + 1);

        if (B0->length - 1 + B1->length - 1 > A->length - 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    _n_fq_neg(alpha, alpha, d, ctx->mod);
    n_fq_bpoly_taylor_shift_gen0_n_fq(B0, alpha, ctx);
    n_fq_bpoly_taylor_shift_gen0_n_fq(B1, alpha, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_fq_set_fq_nmod(alpha, alpha_, ctx);
        _n_fq_neg(alpha, alpha, d, ctx->mod);
        n_fq_bpoly_taylor_shift_gen0_n_fq(A, alpha, ctx);
        n_fq_bpoly_mul(tp1, B0, B1, ctx);
        FLINT_ASSERT(n_fq_bpoly_equal(tp1, A, ctx));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }
#endif

    n_poly_stack_give_back(St->poly_stack, 7);
    n_bpoly_stack_give_back(St->bpoly_stack, 2);

    flint_free(alpha);

    return success;
}


int n_fq_bpoly_hlift2(
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_t B0,
    n_bpoly_t B1,
    const fq_nmod_t alpha_,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx,
    n_poly_bpoly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    int success;
    slong i, j;
    n_fq_poly_struct * c, * s, * t, * u, * v, * g;
    mp_limb_t * alpha;

    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(B0, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(B1, ctx));

    if (A->length < 1 || B0->length < 1 || B1->length < 1)
        return -1;

    n_poly_stack_fit_request(St->poly_stack, 6);
    c = n_poly_stack_take_top(St->poly_stack);
    s = n_poly_stack_take_top(St->poly_stack);
    t = n_poly_stack_take_top(St->poly_stack);
    u = n_poly_stack_take_top(St->poly_stack);
    v = n_poly_stack_take_top(St->poly_stack);
    g = n_poly_stack_take_top(St->poly_stack);

    alpha = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    n_fq_set_fq_nmod(alpha, alpha_, ctx);

    n_fq_bpoly_taylor_shift_gen0_n_fq(A, alpha, ctx);
    n_fq_bpoly_taylor_shift_gen0_n_fq(B0, alpha, ctx);
    n_fq_bpoly_taylor_shift_gen0_n_fq(B1, alpha, ctx);

#if FLINT_WANT_ASSERT
    {
        n_poly_t T;
        n_poly_init(T);
        n_fq_poly_mul(T, B0->coeffs + 0, B1->coeffs + 0, ctx);
        FLINT_ASSERT(n_fq_poly_equal(A->coeffs + 0, T, ctx));
        n_poly_clear(T);
    }
#endif

    if (n_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(n_bpoly_degree1(A) == n_poly_degree(A->coeffs + 0));

    n_fq_poly_xgcd(g, s, t, B1->coeffs + 0, B0->coeffs + 0, ctx);
    if (!n_fq_poly_is_one(g, ctx))
    {
        success = -1;
        goto cleanup;
    }

    n_bpoly_fit_length(B0, A->length);
    n_bpoly_fit_length(B1, A->length);
    for (j = 1; j < A->length; j++)
    {
        n_fq_poly_set(c, A->coeffs + j, ctx);
        for (i = 0; i <= j; i++)
        {
            if (i < B0->length && j - i < B1->length)
            {
                n_fq_poly_mul_(t, B0->coeffs + i, B1->coeffs + j - i, ctx, St->poly_stack);
                n_fq_poly_sub(c, c, t, ctx);
            }
        }

        if (n_fq_poly_is_zero(c))
            continue;

        n_fq_poly_mul_(t, s, c, ctx, St->poly_stack);
        n_fq_poly_divrem_(g, u, t, B0->coeffs + 0, ctx, St->poly_stack);
        n_fq_poly_mul_(t, u, B1->coeffs + 0, ctx, St->poly_stack);
        n_fq_poly_sub(c, c, t, ctx);
        n_fq_poly_divrem_(v, g, c, B0->coeffs + 0, ctx, St->poly_stack);

        if (j < B0->length)
            n_fq_poly_add(B0->coeffs + j, B0->coeffs + j, u, ctx);
        else
            n_fq_poly_set(B0->coeffs + j, u, ctx);

        if (j < B1->length)
            n_fq_poly_add(B1->coeffs + j, B1->coeffs + j, v, ctx);
        else
            n_fq_poly_set(B1->coeffs + j, v, ctx);

        if (!n_fq_poly_is_zero(B0->coeffs + j))
            B0->length = FLINT_MAX(B0->length, j + 1);
        if (!n_fq_poly_is_zero(B1->coeffs + j))
            B1->length = FLINT_MAX(B1->length, j + 1);

        if (B0->length - 1 + B1->length - 1 > A->length - 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    _n_fq_neg(alpha, alpha, d, ctx->mod);
    n_fq_bpoly_taylor_shift_gen0_n_fq(B0, alpha, ctx);
    n_fq_bpoly_taylor_shift_gen0_n_fq(B1, alpha, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_fq_set_fq_nmod(alpha, alpha_, ctx);
        _n_fq_neg(alpha, alpha, d, ctx->mod);
        n_fq_bpoly_taylor_shift_gen0_n_fq(A, alpha, ctx);
        n_fq_bpoly_mul(tp1, B0, B1, ctx);
        FLINT_ASSERT(n_fq_bpoly_equal(tp1, A, ctx));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }
#endif

    n_poly_stack_give_back(St->poly_stack, 6);

    flint_free(alpha);

    return success;
}


/* r factor version */
int n_fq_bpoly_hlift(
    slong r,
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_struct * B,
    const fq_nmod_t alpha,
    slong degree_inner, /* required degree in x */
    const fq_nmod_ctx_t ctx,
    n_poly_bpoly_stack_t St)
{
    int success;
    slong i, j, k, tdeg;
    n_poly_struct * s, * v;
    n_poly_t c, t, u, g1, g2;
    n_bpoly_struct * U;
    fq_nmod_t malpha;

    FLINT_ASSERT(r > 2);

    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, ctx));
    if (A->length < 1)
        return -1;

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(n_fq_bpoly_is_canonical(B + i, ctx));
        if (B[i].length < 1)
            return -1;
    }

    U = FLINT_ARRAY_ALLOC(r, n_bpoly_struct);
    for (i = 0; i < r; i++)
    {
        n_bpoly_init(U + i);
        n_bpoly_fit_length(U + i, A->length);
        for (j = 0; j < A->length; j++)
            n_poly_zero(U[i].coeffs + j);
        U[i].length = A->length;
        n_bpoly_fit_length(B + i, A->length);
    }

    s = FLINT_ARRAY_ALLOC(r, n_poly_struct);
    v = FLINT_ARRAY_ALLOC(r, n_poly_struct);
    for (i = 0; i < r; i++)
    {
        n_poly_init(s + i);
        n_poly_init(v + i);
    }

    n_poly_init(c);
    n_poly_init(t);
    n_poly_init(u);
    n_poly_init(g1);
    n_poly_init(g2);
    fq_nmod_init(malpha, ctx);

    fq_nmod_neg(malpha, alpha, ctx);

    n_fq_bpoly_taylor_shift_gen0_fq_nmod(A, alpha, ctx);
    for (i = 0; i < r; i++)
        n_fq_bpoly_taylor_shift_gen0_fq_nmod(B + i, alpha, ctx);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) * ... */
#if FLINT_WANT_ASSERT
    {
        n_poly_t T;
        n_poly_init(T);
        n_fq_poly_mul(T, B[0].coeffs + 0, B[1].coeffs + 0, ctx);
        for (i = 2; i < r; i++)
            n_fq_poly_mul(T, T, B[i].coeffs + 0, ctx);
        FLINT_ASSERT(n_fq_poly_equal(A->coeffs + 0, T, ctx));
        n_poly_clear(T);
    }
#endif

    if (n_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(n_bpoly_degree1(A) == n_poly_degree(A->coeffs + 0));

    for (i = 0; i < r; i++)
    {
        n_fq_poly_one(t, ctx);
        for (j = 0; j < r; j++)
        {
            if (j != i)
                n_fq_poly_mul(t, t, B[j].coeffs + 0, ctx);
        }

        n_fq_poly_xgcd(g1, s + i, g2, t, B[i].coeffs + 0, ctx);
        if (!n_fq_poly_is_one(g1, ctx))
        {
            success = -1;
            goto cleanup;
        }
    }

    k = r - 2;
    n_fq_poly_mul(U[k].coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, ctx);
    for (k--; k > 0; k--)
        n_fq_poly_mul(U[k].coeffs + 0, B[k].coeffs + 0, U[k + 1].coeffs + 0, ctx);

    for (j = 1; j < A->length; j++)
    {
        for (k = 0; k < r; k++)
            n_poly_zero(U[k].coeffs + j);

        k = r - 2;
        n_poly_zero(U[k].coeffs + j);
        for (i = 0; i <= j; i++)
        {
            if (i < B[k].length && j - i < B[k + 1].length)
            {
                n_fq_poly_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
                n_fq_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
            }
        }
        for (k--; k > 0; k--)
        {
            n_poly_zero(U[k].coeffs + j);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length)
                {
                    n_fq_poly_mul(t, B[k].coeffs + i, U[k + 1].coeffs + j - i, ctx);
                    n_fq_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
                }
            }
        }

        n_fq_poly_set(c, A->coeffs + j, ctx);

        for (i = 0; i <= j; i++)
        {
            if (i < B[0].length)
            {
                n_fq_poly_mul(t, B[0].coeffs + i, U[1].coeffs + j - i, ctx);
                n_fq_poly_sub(c, c, t, ctx);
            }
        }

        if (n_poly_is_zero(c))
            continue;

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            n_fq_poly_mul(t, s + i, c, ctx);
            n_fq_poly_divrem(g1, v + i, t, B[i].coeffs + 0, ctx);
            while (j >= B[i].length)
            {
                n_poly_zero(B[i].coeffs + B[i].length);
                B[i].length++;
            }
            n_fq_poly_add(B[i].coeffs + j, B[i].coeffs + j, v + i, ctx);
            n_bpoly_normalise(B + i);
            tdeg += B[i].length - 1;
        }

        if (tdeg >= A->length)
        {
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        n_fq_poly_mul(t, B[k].coeffs + 0, v + k + 1, ctx);
        n_fq_poly_mul(u, B[k + 1].coeffs + 0, v + k, ctx);
        n_fq_poly_add(t, t, u, ctx);
        n_fq_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
        for (k--; k > 0; k--)
        {
            n_fq_poly_mul(u, B[k].coeffs + 0, t, ctx);
            n_fq_poly_mul(t, U[k + 1].coeffs + 0, v + k, ctx);
            n_fq_poly_add(t, t, u, ctx);
            n_fq_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, ctx);
        }
    }

    for (i = 0; i < r; i++)
        n_fq_bpoly_taylor_shift_gen0_fq_nmod(B + i, malpha, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        n_fq_bpoly_t tp1, tp2;
        n_fq_bpoly_init(tp1);
        n_fq_bpoly_init(tp2);

        n_fq_bpoly_taylor_shift_gen0_fq_nmod(A, malpha, ctx);
        n_fq_bpoly_mul(tp1, B + 0, B + 1, ctx);
        for (i = 2; i < r; i++)
        {
            n_fq_bpoly_mul(tp2, tp1, B + i, ctx);
            n_bpoly_swap(tp1, tp2);
        }
        FLINT_ASSERT(n_fq_bpoly_equal(tp1, A, ctx));

        n_fq_bpoly_clear(tp1);
        n_fq_bpoly_clear(tp2);
    }
#endif

    for (i = 0; i < r; i++)
    {
        n_bpoly_clear(U + i);
        n_poly_clear(s + i);
        n_poly_clear(v + i);
    }
    flint_free(U);
    flint_free(s);
    flint_free(v);

    n_poly_clear(c);
    n_poly_clear(t);
    n_poly_clear(u);
    n_poly_clear(g1);
    n_poly_clear(g2);

    fq_nmod_clear(malpha, ctx);

    return success;
}
