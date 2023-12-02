/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "fmpz_mat.h"
#include "bool_mat.h"
#include "arb_mat.h"

#define LOG2_OVER_E 0.25499459743395350926

slong
_arb_mat_exp_choose_N(const mag_t norm, slong prec)
{
    if (mag_is_special(norm) || mag_cmp_2exp_si(norm, 30) > 0)
    {
        return 1;
    }
    else if (mag_cmp_2exp_si(norm, -prec) < 0)
    {
        return 2;
    }
    else if (mag_cmp_2exp_si(norm, -300) < 0)
    {
        slong N = -MAG_EXP(norm);
        return (prec + N - 1) / N;
    }
    else
    {
        double c, t;

        c = mag_get_d(norm);
        t = d_lambertw(prec * LOG2_OVER_E / c);
        t = c * exp(t + 1.0);
        return FLINT_MIN((slong) (t + 1.0), 2 * prec);
    }
}
static void
_arb_mat_exp_diagonal(arb_mat_t B, const arb_mat_t A, slong prec)
{
    slong n, i;
    n = arb_mat_nrows(A);
    if (B != A)
    {
        arb_mat_zero(B);
    }
    for (i = 0; i < n; i++)
    {
        arb_exp(arb_mat_entry(B, i, i), arb_mat_entry(A, i, i), prec);
    }
}

void
arb_mat_exp(arb_mat_t B, const arb_mat_t A, slong prec)
{
    slong i, j, dim, nz;
    bool_mat_t S;
    slong nildegree;

    if (!arb_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "arb_mat_exp: a square matrix is required!\n");
    }

    if (arb_mat_is_empty(A))
        return;

    dim = arb_mat_nrows(A);

    if (dim == 1)
    {
        arb_exp(arb_mat_entry(B, 0, 0), arb_mat_entry(A, 0, 0), prec);
        return;
    }

    nz = arb_mat_count_is_zero(A);

    if (nz == dim * dim)
    {
        arb_mat_one(B);
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
                bool_mat_set_entry(S, i, j, !arb_is_zero(arb_mat_entry(A, i, j)));
        if (bool_mat_is_diagonal(S))
        {
            _arb_mat_exp_diagonal(B, A, prec);
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
        arb_mat_t T;

        wp = prec + 3 * FLINT_BIT_COUNT(prec);

        mag_init(norm);
        mag_init(err);
        arb_mat_init(T, dim, dim);

        arb_mat_bound_inf_norm(norm, A);

        q = pow(wp, 0.25);  /* wanted magnitude */

        if (mag_cmp_2exp_si(norm, 2 * wp) > 0) /* too big */
            r = 2 * wp;
        else if (mag_cmp_2exp_si(norm, -q) < 0) /* tiny, no need to reduce */
            r = 0;
        else
            r = FLINT_MAX(0, q + MAG_EXP(norm)); /* reduce to magnitude 2^(-r) */

        arb_mat_scalar_mul_2exp_si(T, A, -r);
        mag_mul_2exp_si(norm, norm, -r);

        N = _arb_mat_exp_choose_N(norm, wp);
        if (N < 1)
            flint_throw(FLINT_ERROR, "(%s): N < 1", __func__);

        /* if positive, nildegree is an upper bound on nilpotency degree */
        if (nildegree > 0)
            N = FLINT_MIN(N, nildegree);

        mag_exp_tail(err, norm, N);
        arb_mat_exp_taylor_sum(B, T, N, wp);

        /* add truncation error to entries for which it is not ruled out */
        if (nz == 0)
        {
            for (i = 0; i < dim; i++)
                for (j = 0; j < dim; j++)
                    arb_add_error_mag(arb_mat_entry(B, i, j), err);
        }
        else if (nildegree < 0 || N < nildegree)
        {
            slong w;
            fmpz_mat_t W;
            fmpz_mat_init(W, dim, dim);
            w = bool_mat_all_pairs_longest_walk(W, S);
            if (w + 1 != nildegree)
                flint_throw(FLINT_ERROR, "(%s): w + 1 != nildegree", __func__);
            for (i = 0; i < dim; i++)
            {
                for (j = 0; j < dim; j++)
                {
                    slong d = fmpz_get_si(fmpz_mat_entry(W, i, j)) + 1;
                    if (d < 0 || N < d)
                    {
                        arb_add_error_mag(arb_mat_entry(B, i, j), err);
                    }
                }
            }
            fmpz_mat_clear(W);
        }

        for (i = 0; i < r; i++)
        {
            arb_mat_sqr(T, B, wp);
            arb_mat_swap(T, B);
        }

        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                arb_set_round(arb_mat_entry(B, i, j),
                    arb_mat_entry(B, i, j), prec);

        mag_clear(norm);
        mag_clear(err);
        arb_mat_clear(T);
    }

    bool_mat_clear(S);
}
