/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int n_bpoly_mod_hlift2(
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_t B0,
    n_bpoly_t B1,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t mod)
{
    int success;
    slong i, j;
    n_poly_t c, s, t, u, v;

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B0, mod));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B1, mod));
    if (A->length < 1 || B0->length < 1 || B1->length < 1)
        return -1;

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B0, mod));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B1, mod));

    n_poly_init(c);
    n_poly_init(s);
    n_poly_init(t);
    n_poly_init(u);
    n_poly_init(v);

    n_bpoly_mod_taylor_shift_var0(A, alpha, mod);
    n_bpoly_mod_taylor_shift_var0(B0, alpha, mod);
    n_bpoly_mod_taylor_shift_var0(B1, alpha, mod);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) */
    FLINT_ASSERT(n_poly_degree(A->coeffs + 0) == n_poly_degree(B0->coeffs + 0) +
                                                 n_poly_degree(B1->coeffs + 0));

    if (n_poly_degree(A->coeffs + 0) != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    /* the required degree in x is supposed to be deg_x(A) */
    FLINT_ASSERT(n_bpoly_degree1(A) == n_poly_degree(A->coeffs + 0));

    if (!n_poly_mod_invmod(s, B1->coeffs + 0, B0->coeffs + 0, mod))
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
                n_poly_mod_mul(t, B0->coeffs + i, B1->coeffs + j - i, mod);
                n_poly_mod_sub(c, c, t, mod);
            }
        }

        if (n_poly_is_zero(c))
            continue;

        n_poly_mod_mul(t, s, c, mod);
        n_poly_mod_rem(u, t, B0->coeffs + 0, mod);
        n_poly_mod_mul(t, u, B1->coeffs + 0, mod);
        n_poly_mod_sub(c, c, t, mod);
        n_poly_mod_div(v, c, B0->coeffs + 0, mod);

        if (j < B0->length)
            n_poly_mod_add(B0->coeffs + j, B0->coeffs + j, u, mod);
        else
            n_poly_set(B0->coeffs + j, u);

        if (j < B1->length)
            n_poly_mod_add(B1->coeffs + j, B1->coeffs + j, v, mod);
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

    n_bpoly_mod_taylor_shift_var0(B0, nmod_neg(alpha, mod), mod);
    n_bpoly_mod_taylor_shift_var0(B1, nmod_neg(alpha, mod), mod);

    success = 1;

cleanup:

    if (success > 0)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_bpoly_mod_taylor_shift_var0(A, nmod_neg(alpha, mod), mod);
        n_bpoly_mod_mul(tp1, B0, B1, mod);
        FLINT_ASSERT(n_bpoly_equal(tp1, A));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }

    n_poly_clear(c);
    n_poly_clear(s);
    n_poly_clear(t);
    n_poly_clear(u);
    n_poly_clear(v);

    return success;
}

/*
    input A, B0, B1 with A(y,x) = B0(y,x) * B1(y,x) mod (y-alpha)
    return
       -1: B0(alpha,x) & B1(alpha,x) are not pairwise prime, or
           A(alpha,x) has wrong degree w.r.t x
        0: lift of B0 and B1 to true factors is impossible
        1: successfully lifted B0 and B1 to true factors without changing lc_x
*/
int n_bpoly_mod_hlift(
    slong r,
    n_bpoly_t A, /* clobbered (shifted by alpha) */
    n_bpoly_struct * B,
    mp_limb_t alpha,
    slong degree_inner, /* required degree in x */
    nmod_t mod)
{
    int success;
    slong i, j, k, tdeg;
    n_poly_struct * s, * v;
    n_poly_t c, t, u;
    n_bpoly_struct * U;

    FLINT_ASSERT(r > 2);

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));
    if (A->length < 1)
        return -1;

    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(n_bpoly_mod_is_canonical(B + i, mod));
        if (B[i].length < 1)
            return -1;
    }

    U = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        n_bpoly_init(U + i);
        n_bpoly_fit_length(U + i, A->length);
        for (j = 0; j < A->length; j++)
            n_poly_zero(U[i].coeffs + j);
        U[i].length = A->length;
        n_bpoly_fit_length(B + i, A->length);
    }

    s = (n_poly_struct *) flint_malloc(r*sizeof(n_poly_struct));
    v = (n_poly_struct *) flint_malloc(r*sizeof(n_poly_struct));
    for (i = 0; i < r; i++)
    {
        n_poly_init(s + i);
        n_poly_init(v + i);
    }

    n_poly_init(c);
    n_poly_init(t);
    n_poly_init(u);

    n_bpoly_mod_taylor_shift_var0(A, alpha, mod);
    for (i = 0; i < r; i++)
        n_bpoly_mod_taylor_shift_var0(B + i, alpha, mod);

    /* supposed to have A(alpha,x) = B0(alpha,x) * B1(alpha,x) * ... */
    j = 0;
    for (i = 0; i < r; i++)
        j += n_poly_degree(B[i].coeffs + 0);
    FLINT_ASSERT(j == n_poly_degree(A->coeffs + 0));

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
                n_poly_mod_mul(t, t, B[j].coeffs + 0, mod);
        }
        if (!n_poly_mod_invmod(s + i, t, B[i].coeffs + 0, mod))
        {
            success = -1;
            goto cleanup;
        }
    }

    k = r - 2;
    n_poly_mod_mul(U[k].coeffs + 0, B[k].coeffs + 0, B[k + 1].coeffs + 0, mod);
    for (k--; k > 0; k--)
        n_poly_mod_mul(U[k].coeffs + 0, B[k].coeffs + 0, U[k + 1].coeffs + 0, mod);

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
                n_poly_mod_mul(t, B[k].coeffs + i, B[k + 1].coeffs + j - i, mod);
                n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, mod);
            }
        }
        for (k--; k > 0; k--)
        {
            n_poly_zero(U[k].coeffs + j);
            for (i = 0; i <= j; i++)
            {
                if (i < B[k].length)
                {
                    n_poly_mod_mul(t, B[k].coeffs + i, U[k + 1].coeffs + j - i, mod);
                    n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, mod);
                }
            }
        }

        n_poly_set(c, A->coeffs + j);

        for (i = 0; i <= j; i++)
        {
            if (i < B[0].length)
            {
                n_poly_mod_mul(t, B[0].coeffs + i, U[1].coeffs + j - i, mod);
                n_poly_mod_sub(c, c, t, mod);
            }
        }

        if (n_poly_is_zero(c))
            continue;

        tdeg = 0;
        for (i = 0; i < r; i++)
        {
            n_poly_mod_mul(t, s + i, c, mod);
            n_poly_mod_rem(v + i, t, B[i].coeffs + 0, mod);
            while (j >= B[i].length)
            {
                n_poly_zero(B[i].coeffs + B[i].length);
                B[i].length++;
            }
            n_poly_mod_add(B[i].coeffs + j, B[i].coeffs + j, v + i, mod);
            n_bpoly_normalise(B + i);
            tdeg += B[i].length - 1;
        }

        if (tdeg >= A->length)
        {
            success = 0;
            goto cleanup;
        }

        k = r - 2;
        n_poly_mod_mul(t, B[k].coeffs + 0, v + k + 1, mod);
        n_poly_mod_mul(u, B[k + 1].coeffs + 0, v + k, mod);
        n_poly_mod_add(t, t, u, mod);
        n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, mod);
        for (k--; k > 0; k--)
        {
            n_poly_mod_mul(u, B[k].coeffs + 0, t, mod);
            n_poly_mod_mul(t, U[k + 1].coeffs + 0, v + k, mod);
            n_poly_mod_add(t, t, u, mod);
            n_poly_mod_add(U[k].coeffs + j, U[k].coeffs + j, t, mod);
        }
    }

    for (i = 0; i < r; i++)
        n_bpoly_mod_taylor_shift_var0(B + i, nmod_neg(alpha, mod), mod);

    success = 1;

cleanup:

    if (success > 0)
    {
        n_bpoly_t tp1, tp2;
        n_bpoly_init(tp1);
        n_bpoly_init(tp2);

        n_bpoly_mod_taylor_shift_var0(A, nmod_neg(alpha, mod), mod);
        n_bpoly_mod_mul(tp1, B + 0, B + 1, mod);
        for (i = 2; i < r; i++)
        {
            n_bpoly_mod_mul(tp2, tp1, B + i, mod);
            n_bpoly_swap(tp1, tp2);
        }
        FLINT_ASSERT(n_bpoly_equal(tp1, A));

        n_bpoly_clear(tp1);
        n_bpoly_clear(tp2);
    }

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

    return success;
}
