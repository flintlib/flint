/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"

int fmpz_mod_bpoly_hlift2(
    fmpz_mod_bpoly_t A, /* clobbered (shifted by alpha) */
    fmpz_mod_bpoly_t B0,
    fmpz_mod_bpoly_t B1,
    const fmpz_t alpha,
    slong degree_inner, /* required degree in x */
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_bpoly_stack_t St)
{
    int success;
    slong i, j;
    fmpz_t malpha;
    fmpz_mod_poly_struct * c, * s, * t, * u, * v;

    FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(A, ctx));
    FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(B0, ctx));
    FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(B1, ctx));

    if (A->length < 1 || B0->length < 1 || B1->length < 1)
        return -1;

    fmpz_init(malpha);

    fmpz_mod_poly_stack_fit_request(St->poly_stack, 5);
    c = fmpz_mod_poly_stack_take_top(St->poly_stack);
    s = fmpz_mod_poly_stack_take_top(St->poly_stack);
    t = fmpz_mod_poly_stack_take_top(St->poly_stack);
    u = fmpz_mod_poly_stack_take_top(St->poly_stack);
    v = fmpz_mod_poly_stack_take_top(St->poly_stack);

    fmpz_mod_bpoly_taylor_shift_gen0(A, alpha, ctx);
    fmpz_mod_bpoly_taylor_shift_gen0(B0, alpha, ctx);
    fmpz_mod_bpoly_taylor_shift_gen0(B1, alpha, ctx);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) */
    FLINT_ASSERT(fmpz_mod_poly_degree(A->coeffs + 0, ctx) ==
                           fmpz_mod_poly_degree(B0->coeffs + 0, ctx) +
                           fmpz_mod_poly_degree(B1->coeffs + 0, ctx));

    if (fmpz_mod_poly_degree(A->coeffs + 0, ctx) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(fmpz_mod_bpoly_degree1(A, ctx) ==
                                     fmpz_mod_poly_degree(A->coeffs + 0, ctx));

    if (!fmpz_mod_poly_invmod(s, B1->coeffs + 0, B0->coeffs + 0, ctx))
    {
        success = -2;
        goto cleanup;
    }

    fmpz_mod_bpoly_fit_length(B0, A->length, ctx);
    fmpz_mod_bpoly_fit_length(B1, A->length, ctx);
    for (j = 1; j < A->length; j++)
    {
        fmpz_mod_poly_set(c, A->coeffs + j, ctx);
        for (i = 0; i <= j; i++)
        {
            if (i < B0->length && j - i < B1->length)
            {
                fmpz_mod_poly_mul(t, B0->coeffs + i, B1->coeffs + j - i, ctx);
                fmpz_mod_poly_sub(c, c, t, ctx);
            }
        }

        if (fmpz_mod_poly_is_zero(c, ctx))
            continue;

        fmpz_mod_poly_mul(t, s, c, ctx);
        fmpz_mod_poly_rem(u, t, B0->coeffs + 0, ctx);
        fmpz_mod_poly_mul(t, u, B1->coeffs + 0, ctx);
        fmpz_mod_poly_sub(c, c, t, ctx);
        fmpz_mod_poly_div(v, c, B0->coeffs + 0, ctx);

        if (j < B0->length)
            fmpz_mod_poly_add(B0->coeffs + j, B0->coeffs + j, u, ctx);
        else
            fmpz_mod_poly_set(B0->coeffs + j, u, ctx);

        if (j < B1->length)
            fmpz_mod_poly_add(B1->coeffs + j, B1->coeffs + j, v, ctx);
        else
            fmpz_mod_poly_set(B1->coeffs + j, v, ctx);

        if (!fmpz_mod_poly_is_zero(B0->coeffs + j, ctx))
            B0->length = FLINT_MAX(B0->length, j + 1);
        if (!fmpz_mod_poly_is_zero(B1->coeffs + j, ctx))
            B1->length = FLINT_MAX(B1->length, j + 1);

        if (B0->length - 1 + B1->length - 1 > A->length - 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    fmpz_mod_neg(malpha, alpha, ctx);
    fmpz_mod_bpoly_taylor_shift_gen0(B0, malpha, ctx);
    fmpz_mod_bpoly_taylor_shift_gen0(B1, malpha, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        fmpz_mod_bpoly_t tp1, tp2;
        fmpz_mod_bpoly_init(tp1, ctx);
        fmpz_mod_bpoly_init(tp2, ctx);

        fmpz_mod_neg(malpha, alpha, ctx);
        fmpz_mod_bpoly_taylor_shift_gen0(A, malpha, ctx);
        fmpz_mod_bpoly_mul(tp1, B0, B1, ctx);
        FLINT_ASSERT(fmpz_mod_bpoly_equal(tp1, A, ctx));

        fmpz_mod_bpoly_clear(tp1, ctx);
        fmpz_mod_bpoly_clear(tp2, ctx);
    }
#endif

    fmpz_clear(malpha);
    fmpz_mod_poly_stack_give_back(St->poly_stack, 5);

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
int fmpz_mod_bpoly_hlift(
    slong r,
    fmpz_mod_bpoly_t A, /* clobbered (shifted by alpha) */
    fmpz_mod_bpoly_struct * B,
    const fmpz_t alpha,
    slong degree_inner, /* required degree in x */
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_bpoly_stack_t St)
{
    int success;
    slong i, j, k, tdeg;
    fmpz_t malpha;
    fmpz_mod_poly_struct * c, * t, * u;
    fmpz_mod_poly_struct ** s, ** v, ** Binv;
    fmpz_mod_bpoly_struct ** U;
    TMP_INIT;

    FLINT_ASSERT(r >= 2);

    if (r < 3)
        return fmpz_mod_bpoly_hlift2(A, B + 0, B + 1, alpha, degree_inner, ctx, St);

    FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(A, ctx));
    if (A->length < 1)
        return -1;

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(B + i, ctx));
        if (B[i].length < 1)
            return -1;
    }

    TMP_START;

    fmpz_mod_bpoly_stack_fit_request(St->bpoly_stack, r);
    U = TMP_ARRAY_ALLOC(r, fmpz_mod_bpoly_struct *);
    for (i = 0; i < r; i++)
    {
        U[i] = fmpz_mod_bpoly_stack_take_top(St->bpoly_stack);
        fmpz_mod_bpoly_fit_length(U[i], A->length, ctx);
        for (j = 0; j < A->length; j++)
            fmpz_mod_poly_zero(U[i]->coeffs + j, ctx);
        U[i]->length = A->length;
        fmpz_mod_bpoly_fit_length(B + i, A->length, ctx);
    }

    fmpz_mod_poly_stack_fit_request(St->poly_stack, 3*r + 3);
    s = TMP_ARRAY_ALLOC(3*r, fmpz_mod_poly_struct *);
    v = s + r;
    Binv = v + r;
    for (i = 0; i < r; i++)
    {
        s[i] = fmpz_mod_poly_stack_take_top(St->poly_stack);
        v[i] = fmpz_mod_poly_stack_take_top(St->poly_stack);
        Binv[i] = fmpz_mod_poly_stack_take_top(St->poly_stack);
    }

    fmpz_init(malpha);
    c = fmpz_mod_poly_stack_take_top(St->poly_stack);
    t = fmpz_mod_poly_stack_take_top(St->poly_stack);
    u = fmpz_mod_poly_stack_take_top(St->poly_stack);

    fmpz_mod_bpoly_taylor_shift_gen0(A, alpha, ctx);
    for (i = 0; i < r; i++)
        fmpz_mod_bpoly_taylor_shift_gen0(B + i, alpha, ctx);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) * ... */
    j = 0;
    for (i = 0; i < r; i++)
        j += fmpz_mod_poly_degree(B[i].coeffs + 0, ctx);
    FLINT_ASSERT(j == fmpz_mod_poly_degree(A->coeffs + 0, ctx));

    if (fmpz_mod_poly_degree(A->coeffs + 0, ctx) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(fmpz_mod_bpoly_degree1(A, ctx) ==
                                     fmpz_mod_poly_degree(A->coeffs + 0, ctx));

    for (i = 0; i < r; i++)
    {
        fmpz_mod_poly_one(t, ctx);
        for (j = 0; j < r; j++)
        {
            if (j != i)
                fmpz_mod_poly_mul(t, t, B[j].coeffs + 0, ctx);
        }

        if (!fmpz_mod_poly_invmod(s[i], t, B[i].coeffs + 0, ctx))
        {
            success = -1;
            goto cleanup;
        }

        fmpz_mod_poly_reverse(t, B[i].coeffs + 0, B[i].coeffs[0].length, ctx);
        fmpz_mod_poly_inv_series(Binv[i], t, B[i].coeffs[0].length, ctx);
    }

    k = r - 2;
    fmpz_mod_poly_mul(U[k]->coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, ctx);
    for (k--; k > 0; k--)
        fmpz_mod_poly_mul(U[k]->coeffs + 0, B[k].coeffs + 0, U[k + 1]->coeffs + 0, ctx);

    for (j = 1; j < A->length; j++)
    {
        for (k = 0; k < r; k++)
            fmpz_mod_poly_zero(U[k]->coeffs + j, ctx);

        k = r - 2;
        fmpz_mod_poly_zero(U[k]->coeffs + j, ctx);
        for (i = 0; i <= j; i++)
        {
            if (i < B[k].length && j - i < B[k + 1].length)
            {
                fmpz_mod_poly_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, ctx);
                fmpz_mod_poly_add(U[k]->coeffs + j, U[k]->coeffs + j, t, ctx);
            }
        }
        for (k--; k > 0; k--)
        {
            fmpz_mod_poly_zero(U[k]->coeffs + j, ctx);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length)
                {
                    fmpz_mod_poly_mul(t, B[k].coeffs + i, U[k + 1]->coeffs + j - i, ctx);
                    fmpz_mod_poly_add(U[k]->coeffs + j, U[k]->coeffs + j, t, ctx);
                }
            }
        }

        fmpz_mod_poly_set(c, A->coeffs + j, ctx);

        for (i = 0; i <= j; i++)
        {
            if (i < B[0].length)
            {
                fmpz_mod_poly_mul(t, B[0].coeffs + i, U[1]->coeffs + j - i, ctx);
                fmpz_mod_poly_sub(c, c, t, ctx);
            }
        }

        if (fmpz_mod_poly_is_zero(c, ctx))
            continue;

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            fmpz_mod_poly_rem(t, c, B[i].coeffs + 0, ctx);
            fmpz_mod_poly_mulmod_preinv(v[i], s[i], t, B[i].coeffs + 0, Binv[i], ctx);
            while (j >= B[i].length)
            {
                fmpz_mod_poly_zero(B[i].coeffs + B[i].length, ctx);
                B[i].length++;
            }
            fmpz_mod_poly_add(B[i].coeffs + j, B[i].coeffs + j, v[i], ctx);
            fmpz_mod_bpoly_normalise(B + i, ctx);
            tdeg += B[i].length - 1;
        }

        if (tdeg >= A->length)
        {
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        fmpz_mod_poly_mul(t, B[k].coeffs + 0, v[k + 1], ctx);
        fmpz_mod_poly_mul(u, B[k + 1].coeffs + 0, v[k], ctx);
        fmpz_mod_poly_add(t, t, u, ctx);
        fmpz_mod_poly_add(U[k]->coeffs + j, U[k]->coeffs + j, t, ctx);
        for (k--; k > 0; k--)
        {
            fmpz_mod_poly_mul(u, B[k].coeffs + 0, t, ctx);
            fmpz_mod_poly_mul(t, U[k + 1]->coeffs + 0, v[k], ctx);
            fmpz_mod_poly_add(t, t, u, ctx);
            fmpz_mod_poly_add(U[k]->coeffs + j, U[k]->coeffs + j, t, ctx);
        }
    }

    fmpz_mod_neg(malpha, alpha, ctx);
    for (i = 0; i < r; i++)
        fmpz_mod_bpoly_taylor_shift_gen0(B + i, malpha, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success > 0)
    {
        fmpz_mod_bpoly_t tp1, tp2;
        fmpz_mod_bpoly_init(tp1, ctx);
        fmpz_mod_bpoly_init(tp2, ctx);

        fmpz_mod_neg(malpha, alpha, ctx);
        fmpz_mod_bpoly_taylor_shift_gen0(A, malpha, ctx);
        fmpz_mod_bpoly_mul(tp1, B + 0, B + 1, ctx);
        for (i = 2; i < r; i++)
        {
            fmpz_mod_bpoly_mul(tp2, tp1, B + i, ctx);
            fmpz_mod_bpoly_swap(tp1, tp2, ctx);
        }
        FLINT_ASSERT(fmpz_mod_bpoly_equal(tp1, A, ctx));

        fmpz_mod_bpoly_clear(tp1, ctx);
        fmpz_mod_bpoly_clear(tp2, ctx);
    }
#endif

    fmpz_clear(malpha);
    fmpz_mod_bpoly_stack_give_back(St->bpoly_stack, r);
    fmpz_mod_poly_stack_give_back(St->poly_stack, 3*r + 3);

    TMP_END;

    return success;
}
