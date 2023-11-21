/*
    Copyright (C) 2018 arbguest

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
_arb_mat_solve_c(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong m, n;
    arb_mat_t I, R;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    arb_mat_init(I, n, n);
    arb_mat_init(R, n, n);

    arb_mat_one(I);
    result = arb_mat_approx_solve(R, A, I, prec);
    if (result)
    {
        arb_mat_t RA, RB;

        arb_mat_init(RA, n, n);
        arb_mat_init(RB, n, m);

        arb_mat_mul(RA, R, A, prec);
        arb_mat_mul(RB, R, B, prec);

        result = arb_mat_solve_lu(X, RA, RB, prec);

        arb_mat_clear(RA);
        arb_mat_clear(RB);
    }

    arb_mat_clear(I);
    arb_mat_clear(R);

    return result;
}

int
_arb_mat_solve_d(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong m, n;
    arb_mat_t I, R;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    arb_mat_init(I, n, n);
    arb_mat_init(R, n, n);

    arb_mat_one(I);
    result = arb_mat_approx_solve(R, A, I, prec);
    if (result)
    {
        arb_mat_t RA, RB, E;
        mag_t d;

        arb_mat_init(RA, n, n);
        arb_mat_init(RB, n, m);
        arb_mat_init(E, n, n);
        mag_init(d);

        arb_mat_mul(RA, R, A, prec);
        arb_mat_mul(RB, R, B, prec);
        arb_mat_sub(E, I, RA, prec);
        arb_mat_bound_inf_norm(d, E);

        if (mag_cmp_2exp_si(d, 0) < 0)
        {
            slong i, j;
            mag_t e, err;

            mag_init(e);
            mag_init(err);

            mag_geom_series(d, d, 1);
            arb_mat_set(X, RB);

            for (j = 0; j < m; j++)
            {
                mag_zero(err);
                for (i = 0; i < n; i++)
                {
                    arb_get_mag(e, arb_mat_entry(RB, i, j));
                    mag_max(err, err, e);
                }
                mag_mul(err, err, d);
                for (i = 0; i < n; i++)
                {
                    arb_add_error_mag(arb_mat_entry(X, i, j), err);
                }
            }

            mag_clear(e);
            mag_clear(err);
        }
        else
        {
            result = arb_mat_solve_lu(X, RA, RB, prec);
        }

        arb_mat_clear(RA);
        arb_mat_clear(RB);
        arb_mat_clear(E);
        mag_clear(d);
    }

    arb_mat_clear(I);
    arb_mat_clear(R);

    return result;
}

int
arb_mat_solve_precond(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong n = arb_mat_nrows(A);
    slong m = arb_mat_ncols(B);

    if (m < 0.1 * n + 1)
        return _arb_mat_solve_c(X, A, B, prec);
    else
        return _arb_mat_solve_d(X, A, B, prec);
}
