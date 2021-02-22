/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

int n_bpoly_mod_hlift2_cubic(
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_t B0,
    n_bpoly_t B1,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t ctx,
    nmod_eval_interp_t E,
    n_poly_bpoly_stack_t St)
{
    int success;
    slong len = nmod_eval_interp_eval_length(E);
    slong i, j;
    n_poly_struct * c, * s, * t, * u, * v, * ce;
    n_bpoly_struct * B0e, * B1e;

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, ctx));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B0, ctx));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B1, ctx));

    if (A->length < 1 || B0->length < 1 || B1->length < 1)
        return -1;

    n_poly_stack_fit_request(St->poly_stack, 6);
    c  = n_poly_stack_take_top(St->poly_stack);
    s  = n_poly_stack_take_top(St->poly_stack);
    t  = n_poly_stack_take_top(St->poly_stack);
    u  = n_poly_stack_take_top(St->poly_stack);
    v  = n_poly_stack_take_top(St->poly_stack);
    ce = n_poly_stack_take_top(St->poly_stack);

    n_bpoly_stack_fit_request(St->bpoly_stack, 2);
    B0e  = n_bpoly_stack_take_top(St->bpoly_stack);
    B1e  = n_bpoly_stack_take_top(St->bpoly_stack);

    n_bpoly_mod_taylor_shift_gen0(A, alpha, ctx);
    n_bpoly_mod_taylor_shift_gen0(B0, alpha, ctx);
    n_bpoly_mod_taylor_shift_gen0(B1, alpha, ctx);

    /* check that A(alpha,x) = B0(alpha,x) * B1(alpha,x) */
#if FLINT_WANT_ASSERT
    {
        n_poly_t T;
        n_poly_init(T);
        n_poly_mod_mul(T, B0->coeffs + 0, B1->coeffs + 0, ctx);
        FLINT_ASSERT(n_poly_equal(A->coeffs + 0, T));
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

    if (!n_poly_mod_invmod(s, B1->coeffs + 0, B0->coeffs + 0, ctx))
    {
        success = -2;
        goto cleanup;
    }

    n_bpoly_fit_length(B0, A->length);
    n_bpoly_fit_length(B0e, A->length);
    for (i = 0; i < B0->length; i++)
        nmod_eval_interp_from_coeffs_poly(B0e->coeffs + i, B0->coeffs + i, E, ctx);
    for (i = B0->length; i < A->length; i++)
    {
        n_poly_zero(B0->coeffs + i);
        nmod_evals_zero(B0e->coeffs + i);
    }

    n_bpoly_fit_length(B1, A->length);
    n_bpoly_fit_length(B1e, A->length);
    for (i = 0; i < B1->length; i++)
        nmod_eval_interp_from_coeffs_poly(B1e->coeffs + i, B1->coeffs + i, E, ctx);
    for (i = B1->length; i < A->length; i++)
    {
        n_poly_zero(B1->coeffs + i);
        nmod_evals_zero(B1e->coeffs + i);
    }

    for (j = 1; j < A->length; j++)
    {
        nmod_evals_zero(ce);
        for (i = 0; i <= j; i++)
        {
            if (i < B0->length && j - i < B1->length)
            {
                nmod_evals_addmul(ce, B0e->coeffs + i,
                                      B1e->coeffs + j - i, len, ctx);
            }
        }

        nmod_eval_interp_to_coeffs_poly(c, ce, E, ctx);
        n_poly_mod_sub(c, A->coeffs + j, c, ctx);

#if FLINT_WANT_ASSERT
        {
            n_poly_t c_check;
            n_poly_init(c_check);
            n_poly_set(c_check, A->coeffs + j);
            for (i = 0; i <= j; i++)
            {
                if (i < B0->length && j - i < B1->length)
                {
                    n_poly_mod_mul(t, B0->coeffs + i, B1->coeffs + j - i, ctx);
                    n_poly_mod_sub(c_check, c_check, t, ctx);
                }
            }

            FLINT_ASSERT(n_poly_equal(c, c_check));

            n_poly_clear(c_check);
        }
#endif

        if (n_poly_is_zero(c))
            continue;

        n_poly_mod_mul(t, s, c, ctx);
        n_poly_mod_rem(u, t, B0->coeffs + 0, ctx);
        n_poly_mod_mul(t, u, B1->coeffs + 0, ctx);
        n_poly_mod_sub(c, c, t, ctx);
        n_poly_mod_div(v, c, B0->coeffs + 0, ctx);

        if (!n_poly_is_zero(u))
        {
            n_poly_mod_add(B0->coeffs + j, B0->coeffs + j, u, ctx);
            nmod_eval_interp_from_coeffs_poly(B0e->coeffs + j, B0->coeffs + j, E, ctx);
        }

        if (!n_poly_is_zero(v))
        {
            n_poly_mod_add(B1->coeffs + j, B1->coeffs + j, v, ctx);
            nmod_eval_interp_from_coeffs_poly(B1e->coeffs + j, B1->coeffs + j, E, ctx);
        }

        if (!n_poly_is_zero(B0->coeffs + j))
            B0->length = FLINT_MAX(B0->length, j + 1);

        if (!n_poly_is_zero(B1->coeffs + j))
            B1->length = FLINT_MAX(B1->length, j + 1);

        if (B0->length - 1 + B1->length - 1 > A->length - 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    n_bpoly_mod_taylor_shift_gen0(B0, nmod_neg(alpha, ctx), ctx);
    n_bpoly_mod_taylor_shift_gen0(B1, nmod_neg(alpha, ctx), ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_bpoly_mod_taylor_shift_gen0(A, nmod_neg(alpha, ctx), ctx);
        n_bpoly_mod_mul(tp1, B0, B1, ctx);
        FLINT_ASSERT(n_bpoly_equal(tp1, A));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }
#endif

    n_poly_stack_give_back(St->poly_stack, 6);
    n_bpoly_stack_give_back(St->bpoly_stack, 2);

    return success;
}


int n_bpoly_mod_hlift2(
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_t B0,
    n_bpoly_t B1,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t ctx,
    n_poly_bpoly_stack_t St)
{
    int success;
    slong i, j;
    n_poly_struct * c, * s, * t, * u, * v;

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, ctx));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B0, ctx));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B1, ctx));

    if (A->length < 1 || B0->length < 1 || B1->length < 1)
        return -1;

    n_poly_stack_fit_request(St->poly_stack, 5);
    c = n_poly_stack_take_top(St->poly_stack);
    s = n_poly_stack_take_top(St->poly_stack);
    t = n_poly_stack_take_top(St->poly_stack);
    u = n_poly_stack_take_top(St->poly_stack);
    v = n_poly_stack_take_top(St->poly_stack);

    n_bpoly_mod_taylor_shift_gen0(A, alpha, ctx);
    n_bpoly_mod_taylor_shift_gen0(B0, alpha, ctx);
    n_bpoly_mod_taylor_shift_gen0(B1, alpha, ctx);

    /* check that A(alpha,x) = B0(alpha,x) * B1(alpha,x) */
#if FLINT_WANT_ASSERT
    {
        n_poly_t T;
        n_poly_init(T);
        n_poly_mod_mul(T, B0->coeffs + 0, B1->coeffs + 0, ctx);
        FLINT_ASSERT(n_poly_equal(A->coeffs + 0, T));
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

    if (!n_poly_mod_invmod(s, B1->coeffs + 0, B0->coeffs + 0, ctx))
    {
        success = -2;
        goto cleanup;
    }

    n_bpoly_fit_length(B0, A->length);
    n_bpoly_fit_length(B1, A->length);
    for (j = 1; j < A->length; j++)
    {
        n_poly_set(c, A->coeffs + j);
        for (i = 0; i <= j; i++)
        {
            if (i < B0->length && j - i < B1->length)
            {
                n_poly_mod_mul(t, B0->coeffs + i, B1->coeffs + j - i, ctx);
                n_poly_mod_sub(c, c, t, ctx);
            }
        }

        if (n_poly_is_zero(c))
            continue;

        n_poly_mod_mul(t, s, c, ctx);
        n_poly_mod_rem(u, t, B0->coeffs + 0, ctx);
        n_poly_mod_mul(t, u, B1->coeffs + 0, ctx);
        n_poly_mod_sub(c, c, t, ctx);
        n_poly_mod_div(v, c, B0->coeffs + 0, ctx);

        if (j < B0->length)
            n_poly_mod_add(B0->coeffs + j, B0->coeffs + j, u, ctx);
        else
            n_poly_set(B0->coeffs + j, u);

        if (j < B1->length)
            n_poly_mod_add(B1->coeffs + j, B1->coeffs + j, v, ctx);
        else
            n_poly_set(B1->coeffs + j, v);

        if (!n_poly_is_zero(B0->coeffs + j))
            B0->length = FLINT_MAX(B0->length, j + 1);
        if (!n_poly_is_zero(B1->coeffs + j))
            B1->length = FLINT_MAX(B1->length, j + 1);

        if (B0->length - 1 + B1->length - 1 > A->length - 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    n_bpoly_mod_taylor_shift_gen0(B0, nmod_neg(alpha, ctx), ctx);
    n_bpoly_mod_taylor_shift_gen0(B1, nmod_neg(alpha, ctx), ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_bpoly_mod_taylor_shift_gen0(A, nmod_neg(alpha, ctx), ctx);
        n_bpoly_mod_mul(tp1, B0, B1, ctx);
        FLINT_ASSERT(n_bpoly_equal(tp1, A));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }
#endif

    n_poly_stack_give_back(St->poly_stack, 5);

    return success;
}


int n_bpoly_mod_hlift_cubic(
    slong r,
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_struct * B,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t ctx,
    nmod_eval_interp_t E,
    n_poly_bpoly_stack_t St)
{
    int success;
    slong len = nmod_eval_interp_eval_length(E);
    slong i, j, k, tdeg;
    n_poly_struct * p, * c, * t, * ce;
    n_poly_struct ** s, * vk, ** Binv, * vek;
    n_bpoly_struct ** Ue, ** Be;
#if FLINT_WANT_ASSERT
    n_bpoly_t Acopy;
    n_bpoly_struct * Bcqt;
    int successcqt;
#endif
    TMP_INIT;

    FLINT_ASSERT(r > 2);

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, ctx));
    if (A->length < 1)
        return -1;

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(n_bpoly_mod_is_canonical(B + i, ctx));
        if (B[i].length < 1)
            return -1;
    }

    TMP_START;

#if FLINT_WANT_ASSERT
    n_bpoly_init(Acopy);
    n_bpoly_set(Acopy, A);
    Bcqt = FLINT_ARRAY_ALLOC(r, n_bpoly_struct);
    for (i = 0; i < r; i++)
    {
        n_bpoly_init(Bcqt + i);
        n_bpoly_set(Bcqt + i, B + i);
    }
    successcqt = n_bpoly_mod_hlift(r, Acopy, Bcqt, alpha, degree_inner, ctx, St);
#endif

    n_bpoly_stack_fit_request(St->bpoly_stack, 2*r);
    Ue = TMP_ARRAY_ALLOC(2*r, n_bpoly_struct *);
    Be = Ue + r;
    for (i = 0; i < r; i++)
    {
        Ue[i] = n_bpoly_stack_take_top(St->bpoly_stack);
        Be[i] = n_bpoly_stack_take_top(St->bpoly_stack);
    }

    n_poly_stack_fit_request(St->poly_stack, 2*r + 5);
    s = TMP_ARRAY_ALLOC(2*r, n_poly_struct *);
    Binv = s + r;
    for (i = 0; i < r; i++)
    {
        s[i] = n_poly_stack_take_top(St->poly_stack);
        Binv[i] = n_poly_stack_take_top(St->poly_stack);
    }
    vk = n_poly_stack_take_top(St->poly_stack);
    vek = n_poly_stack_take_top(St->poly_stack);
    c = n_poly_stack_take_top(St->poly_stack);
    t = n_poly_stack_take_top(St->poly_stack);
    ce = n_poly_stack_take_top(St->poly_stack);

    n_bpoly_mod_taylor_shift_gen0(A, alpha, ctx);
    for (i = 0; i < r; i++)
        n_bpoly_mod_taylor_shift_gen0(B + i, alpha, ctx);

    /* check that A(alpha,x) = B0(alpha,x) * B1(alpha,x) * ... */
#if FLINT_WANT_ASSERT
    {
        n_poly_t T;
        n_poly_init(T);
        n_poly_mod_mul(T, B[0].coeffs + 0, B[1].coeffs + 0, ctx);
        for (i = 2; i < r; i++)
            n_poly_mod_mul(T, T, B[i].coeffs + 0, ctx);
        FLINT_ASSERT(n_poly_equal(A->coeffs + 0, T));
        n_poly_clear(T);
    }
#endif

    /* the required degree in x is supposed to be deg_x(A) */
    if (n_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    FLINT_ASSERT(n_bpoly_degree1(A) == degree_inner);

    for (k = 0; k < r; k++)
    {
        /* s[k] = (prod_{i!=k} B[i].coeffs[0])^-1 (mod B[k].coeffs[0]) */
        n_poly_mod_div(t, A->coeffs + 0, B[k].coeffs + 0, ctx);
        if (!n_poly_mod_invmod(s[k], t, B[k].coeffs + 0, ctx))
        {
            success = -1;
            goto cleanup;
        }

        /* set up mul (mod B[k].coeffs[0]) */
        n_poly_reverse(t, B[k].coeffs + 0, B[k].coeffs[0].length);
        n_poly_mod_inv_series(Binv[k], t, B[k].coeffs[0].length, ctx);

        /* set up the evaluation of B[k] */
        n_bpoly_fit_length(B + k, A->length);
        n_bpoly_fit_length(Be[k], A->length);

        for (i = 0; i < B[k].length; i++)
            nmod_eval_interp_from_coeffs_poly(Be[k]->coeffs + i,
                                              B[k].coeffs + i, E, ctx);

        for (i = B[k].length; i < A->length; i++)
        {
            n_poly_zero(B[k].coeffs + i);
            nmod_evals_zero(Be[k]->coeffs + i);
        }

        /* Ue[0] is not used */
        if (k > 0)
        {
            n_bpoly_fit_length(Ue[k], A->length);
            Ue[k]->length = A->length;
            for (i = 0; i < A->length; i++)
                nmod_evals_zero(Ue[k]->coeffs + i);
        }
    }

    k = r - 2;
    nmod_evals_mul(Ue[k]->coeffs + 0, Be[k]->coeffs + 0,
                                      Be[k + 1]->coeffs + 0, len, ctx);
    for (k--; k > 0; k--)
        nmod_evals_mul(Ue[k]->coeffs + 0, Be[k]->coeffs + 0,
                                          Ue[k + 1]->coeffs + 0, len, ctx);

    for (j = 1; j < A->length; j++)
    {
        k = r - 2;
        nmod_evals_zero(Ue[k]->coeffs + j);
        for (i = 0; i <= j; i++)
            nmod_evals_addmul(Ue[k]->coeffs + j, Be[k]->coeffs + i,
                                          Be[k + 1]->coeffs + j - i, len, ctx);
        for (k--; k > 0; k--)
        {
            nmod_evals_zero(Ue[k]->coeffs + j);
            for (i = 0; i <= j; i++)
                nmod_evals_addmul(Ue[k]->coeffs + j, Be[k]->coeffs + i,
                                          Ue[k + 1]->coeffs + j - i, len, ctx);
        }

        nmod_evals_zero(ce);
        for (i = 0; i <= j; i++)
            nmod_evals_addmul(ce, Be[0]->coeffs + i, Ue[1]->coeffs + j - i, len, ctx);

        nmod_eval_interp_to_coeffs_poly(c, ce, E, ctx);
        n_poly_mod_sub(c, A->coeffs + j, c, ctx);

        if (n_poly_is_zero(c))
            continue;

        tdeg = 0;
        for (k = r - 1; k >= 0; k--)
        {
            n_poly_mod_rem(t, c, B[k].coeffs + 0, ctx);
            n_poly_mod_mulmod_preinv(vk, s[k], t, B[k].coeffs + 0, Binv[k], ctx);
            nmod_eval_interp_from_coeffs_poly(vek, vk, E, ctx);
            if (!n_poly_is_zero(vk))
            {
                nmod_evals_add_inplace(Be[k]->coeffs + j, vek, len, ctx);
                n_poly_mod_add(B[k].coeffs + j, B[k].coeffs + j, vk, ctx);
                if (!n_poly_is_zero(B[k].coeffs + j))
                    B[k].length = FLINT_MAX(B[k].length, j + 1);
            }
            tdeg += B[k].length - 1;

            /* correct the U's */
            if (k > r - 2)
            {
                n_poly_swap(ce, vek);
            }
            else if (k > 0)
            {
                p = (k == r - 2) ? Be[k + 1]->coeffs : Ue[k + 1]->coeffs;
                nmod_evals_fmma(ce, Be[k]->coeffs + 0, ce, p, vek, len, ctx);
                nmod_evals_add_inplace(Ue[k]->coeffs + j, ce, len, ctx);
            }
        }

        if (tdeg >= A->length)
        {
            success = 0;
            goto cleanup;
        }
    }

    for (k = 0; k < r; k++)
        n_bpoly_mod_taylor_shift_gen0(B + k, nmod_neg(alpha, ctx), ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    FLINT_ASSERT(success == successcqt);
    if (success > 0)
    {
        for (i = 0; i < r; i++)
            FLINT_ASSERT(n_bpoly_equal(Bcqt + i, B + i));
    }

    n_bpoly_clear(Acopy);
    for (i = 0; i < r; i++)
        n_bpoly_clear(Bcqt + i);
    flint_free(Bcqt);
#endif

    n_bpoly_stack_give_back(St->bpoly_stack, 2*r);
    n_poly_stack_give_back(St->poly_stack, 2*r + 5);

    TMP_END;

    return success;
}


/*
    y = gen(0), x = gen(1)
    input A, Bi with A(y,x) = prod_i Bi(y,x) ctx (y-alpha)
    return
       -1: the Bi(alpha,x) are not pairwise prime, or
           A(alpha,x) has wrong degree w.r.t x
        0: lift of the Bi to true factors is impossible
        1: successfully lifted the Bi to true factors without changing lc_x
*/
int n_bpoly_mod_hlift(
    slong r,
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_struct * B,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t ctx,
    n_poly_bpoly_stack_t St)
{
    int success;
    slong i, j, k, tdeg;
    n_poly_struct * c, * t, * u;
    n_poly_struct ** s, ** v, ** Binv;
    n_bpoly_struct ** U;
    TMP_INIT;

    FLINT_ASSERT(r > 2);

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, ctx));
    if (A->length < 1)
        return -1;

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(n_bpoly_mod_is_canonical(B + i, ctx));
        if (B[i].length < 1)
            return -1;
    }

    TMP_START;

    n_bpoly_stack_fit_request(St->bpoly_stack, r);
    U = TMP_ARRAY_ALLOC(r, n_bpoly_struct *);
    for (i = 0; i < r; i++)
    {
        U[i] = n_bpoly_stack_take_top(St->bpoly_stack);
        n_bpoly_fit_length(U[i], A->length);
        for (j = 0; j < A->length; j++)
            n_poly_zero(U[i]->coeffs + j);
        U[i]->length = A->length;
        n_bpoly_fit_length(B + i, A->length);
    }

    n_poly_stack_fit_request(St->poly_stack, 3*r + 3);
    s = TMP_ARRAY_ALLOC(3*r, n_poly_struct *);
    v = s + r;
    Binv = v + r;
    for (i = 0; i < r; i++)
    {
        s[i] = n_poly_stack_take_top(St->poly_stack);
        v[i] = n_poly_stack_take_top(St->poly_stack);
        Binv[i] = n_poly_stack_take_top(St->poly_stack);
    }

    c = n_poly_stack_take_top(St->poly_stack);
    t = n_poly_stack_take_top(St->poly_stack);
    u = n_poly_stack_take_top(St->poly_stack);

    n_bpoly_mod_taylor_shift_gen0(A, alpha, ctx);
    for (i = 0; i < r; i++)
        n_bpoly_mod_taylor_shift_gen0(B + i, alpha, ctx);

    /* check that A(alpha,x) = B0(alpha,x) * B1(alpha,x) * ... */
#if FLINT_WANT_ASSERT
    {
        n_poly_t T;
        n_poly_init(T);
        n_poly_mod_mul(T, B[0].coeffs + 0, B[1].coeffs + 0, ctx);
        for (i = 2; i < r; i++)
            n_poly_mod_mul(T, T, B[i].coeffs + 0, ctx);
        FLINT_ASSERT(n_poly_equal(A->coeffs + 0, T));
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
        n_poly_one(t);
        for (j = 0; j < r; j++)
        {
            if (j != i)
                n_poly_mod_mul(t, t, B[j].coeffs + 0, ctx);
        }

        if (!n_poly_mod_invmod(s[i], t, B[i].coeffs + 0, ctx))
        {
            success = -1;
            goto cleanup;
        }

        n_poly_reverse(t, B[i].coeffs + 0, B[i].coeffs[0].length);
        n_poly_mod_inv_series(Binv[i], t, B[i].coeffs[0].length, ctx);
    }

    k = r - 2;
    n_poly_mod_mul(U[k]->coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, ctx);
    for (k--; k > 0; k--)
        n_poly_mod_mul(U[k]->coeffs + 0, B[k].coeffs + 0, U[k + 1]->coeffs + 0, ctx);

    for (j = 1; j < A->length; j++)
    {
        for (k = 0; k < r; k++)
            n_poly_zero(U[k]->coeffs + j);

        k = r - 2;
        n_poly_zero(U[k]->coeffs + j);
        for (i = 0; i <= j; i++)
        {
            if (i < B[k].length && j - i < B[k + 1].length)
            {
                n_poly_mod_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
                n_poly_mod_add(U[k]->coeffs + j, U[k]->coeffs + j, t, ctx);
            }
        }
        for (k--; k > 0; k--)
        {
            n_poly_zero(U[k]->coeffs + j);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length)
                {
                    n_poly_mod_mul(t, B[k].coeffs + i, U[k + 1]->coeffs + j - i, ctx);
                    n_poly_mod_add(U[k]->coeffs + j, U[k]->coeffs + j, t, ctx);
                }
            }
        }

        n_poly_set(c, A->coeffs + j);

        for (i = 0; i <= j; i++)
        {
            if (i < B[0].length)
            {
                n_poly_mod_mul(t, B[0].coeffs + i, U[1]->coeffs + j - i, ctx);
                n_poly_mod_sub(c, c, t, ctx);
            }
        }

        if (n_poly_is_zero(c))
            continue;

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            n_poly_mod_rem(t, c, B[i].coeffs + 0, ctx);
            n_poly_mod_mulmod_preinv(v[i], s[i], t, B[i].coeffs + 0, Binv[i], ctx);
            while (j >= B[i].length)
            {
                n_poly_zero(B[i].coeffs + B[i].length);
                B[i].length++;
            }
            n_poly_mod_add(B[i].coeffs + j, B[i].coeffs + j, v[i], ctx);
            n_bpoly_normalise(B + i);
            tdeg += B[i].length - 1;
        }

        if (tdeg >= A->length)
        {
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        n_poly_mod_mul(t, B[k].coeffs + 0, v[k + 1], ctx);
        n_poly_mod_mul(u, B[k + 1].coeffs + 0, v[k], ctx);
        n_poly_mod_add(t, t, u, ctx);
        n_poly_mod_add(U[k]->coeffs + j, U[k]->coeffs + j, t, ctx);
        for (k--; k > 0; k--)
        {
            n_poly_mod_mul(u, B[k].coeffs + 0, t, ctx);
            n_poly_mod_mul(t, U[k + 1]->coeffs + 0, v[k], ctx);
            n_poly_mod_add(t, t, u, ctx);
            n_poly_mod_add(U[k]->coeffs + j, U[k]->coeffs + j, t, ctx);
        }
    }

    for (i = 0; i < r; i++)
        n_bpoly_mod_taylor_shift_gen0(B + i, nmod_neg(alpha, ctx), ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_bpoly_mod_taylor_shift_gen0(A, nmod_neg(alpha, ctx), ctx);
        n_bpoly_mod_mul(tp1, B + 0, B + 1, ctx);
        for (i = 2; i < r; i++)
        {
            n_bpoly_mod_mul(tp2, tp1, B + i, ctx);
            n_bpoly_swap(tp1, tp2);
        }
        FLINT_ASSERT(n_bpoly_equal(tp1, A));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }
#endif

    n_bpoly_stack_give_back(St->bpoly_stack, r);
    n_poly_stack_give_back(St->poly_stack, 3*r + 3);

    TMP_END;

    return success;
}
