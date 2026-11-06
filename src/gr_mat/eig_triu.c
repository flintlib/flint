/*
    Copyright (C) 2013 Timo Hartmann
    Copyright (C) 2018, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

/*
    Right eigenvectors of an upper triangular matrix, by back
    substitution with rescaling to avoid overflow (compare LAPACK's
    ztrevc).  Aliasing between ER and A is not allowed.

    The rescaling constants are chosen exactly as in the acb version:
    smlnum = n * 2^(-29*prec) guards the diagonal against near-zero
    divisors, and a running row maximum is divided out whenever it
    exceeds 2^(prec/2).  In a ring with a bounded exponent range (such
    as the double expansions of dfloat) smlnum underflows to zero; it is
    then n * 2^(-prec) times the largest entry of A instead (or n *
    2^(-prec) for a zero matrix): a perturbation of the diagonal at the
    level of the backward error, small enough for accuracy and large
    enough that the quotients stay far from overflow between two
    rescalings.
*/

static int
_gr_triu_constants(gr_ptr smlnum, gr_ptr simin, const gr_mat_t A, slong n, slong prec, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_one(smlnum, ctx);
    status |= gr_mul_2exp_si(smlnum, smlnum, -29 * prec, ctx);
    status |= gr_mul_ui(smlnum, smlnum, n, ctx);

    if (status == GR_SUCCESS && gr_is_zero(smlnum, ctx) == T_TRUE)
    {
        status |= gr_mat_norm_max(smlnum, A, ctx);
        if (status == GR_SUCCESS && gr_is_zero(smlnum, ctx) == T_TRUE)
            status |= gr_one(smlnum, ctx);
        status |= gr_mul_2exp_si(smlnum, smlnum, -prec, ctx);
        status |= gr_mul_ui(smlnum, smlnum, n, ctx);
    }

    status |= gr_one(simin, ctx);
    status |= gr_mul_2exp_si(simin, simin, prec / 2, ctx);

    return status;
}

int
gr_mat_eig_triu_r(gr_mat_t ER, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j, k, n, prec;
    slong sz = ctx->sizeof_elem;
    gr_ptr r, s, t, tm, smin, smlnum, simin, rmax;
    /* tm is retained for the smin computation */
    int status = GR_SUCCESS;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx) || gr_mat_nrows(ER, ctx) != n
            || gr_mat_ncols(ER, ctx) != n)
        return GR_DOMAIN;

    if (gr_ctx_get_real_prec(&prec, ctx) != GR_SUCCESS)
        return GR_UNABLE;

    status |= gr_mat_one(ER, ctx);

    if (n <= 1)
        return status;

    GR_TMP_INIT5(r, s, t, tm, smin, ctx);
    GR_TMP_INIT3(smlnum, simin, rmax, ctx);

    status |= _gr_triu_constants(smlnum, simin, A, n, prec, ctx);

    for (i = 1; i < n; i++)
    {
        status |= gr_one(rmax, ctx);
        status |= gr_set(s, GR_MAT_ENTRY(A, i, i, sz), ctx);

        /* smin = max(eps * |s|, smlnum) */
        status |= gr_abs(smin, s, ctx);
        status |= gr_mul_2exp_si(smin, smin, -prec, ctx);
        if (gr_abs_lt(smin, smlnum, ctx) == T_TRUE)
            status |= gr_set(smin, smlnum, ctx);

        for (j = i - 1; j >= 0; j--)
        {
            status |= _gr_vec_dot(r, NULL, 0, GR_MAT_ENTRY(A, j, j + 1, sz),
                                              GR_MAT_ENTRY(ER, i, j + 1, sz), i - j, ctx);

            status |= gr_sub(t, GR_MAT_ENTRY(A, j, j, sz), s, ctx);

            if (gr_abs_lt(t, smin, ctx) == T_TRUE)
                status |= gr_set(t, smin, ctx);

            status |= gr_div(GR_MAT_ENTRY(ER, i, j, sz), r, t, ctx);
            status |= gr_neg(GR_MAT_ENTRY(ER, i, j, sz), GR_MAT_ENTRY(ER, i, j, sz), ctx);

            /* Track the running maximum as a ring element rather than a
               magnitude: only |rmax| matters for the overflow test, and
               dividing a row by a complex number of the right magnitude
               rescales the eigenvector by a phase, which is immaterial.
               This keeps a square root out of an O(n^2) loop. */
            if (gr_abs_gt(r, rmax, ctx) == T_TRUE)
                status |= gr_set(rmax, r, ctx);

            if (gr_abs_gt(rmax, simin, ctx) == T_TRUE)
            {
                for (k = j; k < i + 1; k++)
                    status |= gr_div(GR_MAT_ENTRY(ER, i, k, sz),
                                     GR_MAT_ENTRY(ER, i, k, sz), rmax, ctx);

                status |= gr_one(rmax, ctx);
            }

            if (status != GR_SUCCESS)
                goto cleanup;
        }

        if (gr_is_one(rmax, ctx) != T_TRUE)
        {
            for (k = 0; k < i + 1; k++)
                status |= gr_div(GR_MAT_ENTRY(ER, i, k, sz),
                                 GR_MAT_ENTRY(ER, i, k, sz), rmax, ctx);
        }
    }

    status |= gr_mat_transpose(ER, ER, ctx);

cleanup:
    GR_TMP_CLEAR5(r, s, t, tm, smin, ctx);
    GR_TMP_CLEAR3(smlnum, simin, rmax, ctx);

    return status;
}

/* Left eigenvectors of an upper triangular matrix.  Aliasing between EL
   and A is not allowed. */
int
gr_mat_eig_triu_l(gr_mat_t EL, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j, k, n, prec;
    slong sz = ctx->sizeof_elem;
    gr_ptr r, s, t, tm, smin, smlnum, simin, rmax;
    /* tm is retained for the smin computation */
    gr_mat_t AT;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx) || gr_mat_nrows(EL, ctx) != n
            || gr_mat_ncols(EL, ctx) != n)
        return GR_DOMAIN;

    if (gr_ctx_get_real_prec(&prec, ctx) != GR_SUCCESS)
        return GR_UNABLE;

    status |= gr_mat_one(EL, ctx);

    if (n <= 1)
        return status;

    gr_mat_init(AT, n, n, ctx);
    status |= gr_mat_transpose(AT, A, ctx);

    GR_TMP_INIT5(r, s, t, tm, smin, ctx);
    GR_TMP_INIT3(smlnum, simin, rmax, ctx);

    status |= _gr_triu_constants(smlnum, simin, A, n, prec, ctx);

    for (i = 0; i < n - 1; i++)
    {
        status |= gr_one(rmax, ctx);
        status |= gr_set(s, GR_MAT_ENTRY(AT, i, i, sz), ctx);

        status |= gr_abs(smin, s, ctx);
        status |= gr_mul_2exp_si(smin, smin, -prec, ctx);
        if (gr_abs_lt(smin, smlnum, ctx) == T_TRUE)
            status |= gr_set(smin, smlnum, ctx);

        for (j = i + 1; j < n; j++)
        {
            status |= _gr_vec_dot(r, NULL, 0, GR_MAT_ENTRY(EL, i, i, sz),
                                              GR_MAT_ENTRY(AT, j, i, sz), j - i, ctx);

            status |= gr_sub(t, GR_MAT_ENTRY(AT, j, j, sz), s, ctx);

            if (gr_abs_lt(t, smin, ctx) == T_TRUE)
                status |= gr_set(t, smin, ctx);

            status |= gr_div(GR_MAT_ENTRY(EL, i, j, sz), r, t, ctx);
            status |= gr_neg(GR_MAT_ENTRY(EL, i, j, sz), GR_MAT_ENTRY(EL, i, j, sz), ctx);

            /* Track the running maximum as a ring element rather than a
               magnitude: only |rmax| matters for the overflow test, and
               dividing a row by a complex number of the right magnitude
               rescales the eigenvector by a phase, which is immaterial.
               This keeps a square root out of an O(n^2) loop. */
            if (gr_abs_gt(r, rmax, ctx) == T_TRUE)
                status |= gr_set(rmax, r, ctx);

            if (gr_abs_gt(rmax, simin, ctx) == T_TRUE)
            {
                for (k = i; k < j + 1; k++)
                    status |= gr_div(GR_MAT_ENTRY(EL, i, k, sz),
                                     GR_MAT_ENTRY(EL, i, k, sz), rmax, ctx);

                status |= gr_one(rmax, ctx);
            }

            if (status != GR_SUCCESS)
                goto cleanup;
        }

        if (gr_is_one(rmax, ctx) != T_TRUE)
        {
            for (k = i; k < n; k++)
                status |= gr_div(GR_MAT_ENTRY(EL, i, k, sz),
                                 GR_MAT_ENTRY(EL, i, k, sz), rmax, ctx);
        }
    }

cleanup:
    gr_mat_clear(AT, ctx);
    GR_TMP_CLEAR5(r, s, t, tm, smin, ctx);
    GR_TMP_CLEAR3(smlnum, simin, rmax, ctx);

    return status;
}
