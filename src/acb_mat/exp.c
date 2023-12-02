/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "fmpz_mat.h"
#include "acb_mat.h"
#include "bool_mat.h"

slong _arb_mat_exp_choose_N(const mag_t norm, slong prec);

static slong
_acb_mat_count_is_zero(const acb_mat_t A)
{
    slong nz, i, j;
    for (nz = 0, i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            nz += acb_is_zero(acb_mat_entry(A, i, j));
    return nz;
}

static void
_acb_mat_exp_diagonal(acb_mat_t B, const acb_mat_t A, slong prec)
{
    slong n, i;
    n = acb_mat_nrows(A);
    if (B != A)
    {
        acb_mat_zero(B);
    }
    for (i = 0; i < n; i++)
    {
        acb_exp(acb_mat_entry(B, i, i), acb_mat_entry(A, i, i), prec);
    }
}

void
acb_mat_exp(acb_mat_t B, const acb_mat_t A, slong prec)
{
    slong i, j, dim, nz;
    bool_mat_t S;
    slong nildegree;

    if (!acb_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "acb_mat_exp: a square matrix is required!\n");
    }

    if (acb_mat_is_empty(A))
        return;

    dim = acb_mat_nrows(A);

    if (dim == 1)
    {
        acb_exp(acb_mat_entry(B, 0, 0), acb_mat_entry(A, 0, 0), prec);
        return;
    }

    if (acb_mat_is_real(A))
    {
        arb_mat_t R;
        arb_mat_init(R, dim, dim);
        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                arb_set(arb_mat_entry(R, i, j),
                        acb_realref(acb_mat_entry(A, i, j)));
        arb_mat_exp(R, R, prec);
        acb_mat_set_arb_mat(B, R);
        arb_mat_clear(R);
        return;
    }

    nz = _acb_mat_count_is_zero(A);

    if (nz == dim * dim)
    {
        acb_mat_one(B);
        return;
    }

    bool_mat_init(S, dim, dim);
    if (nz == 0)
    {
        nildegree = -1;
        bool_mat_complement(S, S);
    }
    else
    {
        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                bool_mat_set_entry(S, i, j, !acb_is_zero(acb_mat_entry(A, i, j)));
        if (bool_mat_is_diagonal(S))
        {
            _acb_mat_exp_diagonal(B, A, prec);
            bool_mat_clear(S);
            return;
        }
        else
        {
            nildegree = bool_mat_nilpotency_degree(S);
        }
    }

    /* evaluate using scaling and squaring of truncated taylor series */
    {
        slong wp, N, q, r;
        mag_t norm, err;
        acb_mat_t T;

        wp = prec + 3 * FLINT_BIT_COUNT(prec);

        mag_init(norm);
        mag_init(err);
        acb_mat_init(T, dim, dim);

        acb_mat_bound_inf_norm(norm, A);

        q = pow(wp, 0.25);  /* wanted magnitude */

        if (mag_cmp_2exp_si(norm, 2 * wp) > 0) /* too big */
            r = 2 * wp;
        else if (mag_cmp_2exp_si(norm, -q) < 0) /* tiny, no need to reduce */
            r = 0;
        else
            r = FLINT_MAX(0, q + MAG_EXP(norm)); /* reduce to magnitude 2^(-r) */

        acb_mat_scalar_mul_2exp_si(T, A, -r);
        mag_mul_2exp_si(norm, norm, -r);

        N = _arb_mat_exp_choose_N(norm, wp);

        /* if positive, nildegree is an upper bound on nilpotency degree */
        if (nildegree > 0)
            N = FLINT_MIN(N, nildegree);

        mag_exp_tail(err, norm, N);
        acb_mat_exp_taylor_sum(B, T, N, wp);

        /* add truncation error to entries for which it is not ruled out */
        if (nz == 0)
        {
            for (i = 0; i < dim; i++)
                for (j = 0; j < dim; j++)
                    acb_add_error_mag(acb_mat_entry(B, i, j), err);
        }
        else if (nildegree < 0 || N < nildegree)
        {
            slong w;
            fmpz_mat_t W;
            fmpz_mat_init(W, dim, dim);
            w = bool_mat_all_pairs_longest_walk(W, S);
            if (w + 1 != nildegree)
                flint_throw(FLINT_ERROR, "(%s)\n", __func__);
            for (i = 0; i < dim; i++)
            {
                for (j = 0; j < dim; j++)
                {
                    slong d = fmpz_get_si(fmpz_mat_entry(W, i, j)) + 1;
                    if (d < 0 || N < d)
                    {
                        acb_add_error_mag(acb_mat_entry(B, i, j), err);
                    }
                }
            }
            fmpz_mat_clear(W);
        }

        for (i = 0; i < r; i++)
        {
            acb_mat_sqr(T, B, wp);
            acb_mat_swap(T, B);
        }

        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                acb_set_round(acb_mat_entry(B, i, j),
                    acb_mat_entry(B, i, j), prec);

        mag_clear(norm);
        mag_clear(err);
        acb_mat_clear(T);
    }

    bool_mat_clear(S);
}
