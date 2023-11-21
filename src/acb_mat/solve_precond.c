/*
    Copyright (C) 2018 arbguest

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
_acb_mat_solve_c(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    int result;
    slong m, n;
    acb_mat_t I, R;

    n = acb_mat_nrows(A);
    m = acb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    acb_mat_init(I, n, n);
    acb_mat_init(R, n, n);

    acb_mat_one(I);
    result = acb_mat_approx_solve(R, A, I, prec);

    if (result)
    {
        acb_mat_t RA, RB;

        acb_mat_init(RA, n, n);
        acb_mat_init(RB, n, m);

        acb_mat_mul(RA, R, A, prec);
        acb_mat_mul(RB, R, B, prec);

        result = acb_mat_solve_lu(X, RA, RB, prec);

        acb_mat_clear(RA);
        acb_mat_clear(RB);
    }

    acb_mat_clear(I);
    acb_mat_clear(R);

    return result;
}

int
_acb_mat_solve_d(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    int result, real;
    slong m, n;
    acb_mat_t I, R;

    n = acb_mat_nrows(A);
    m = acb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    real = acb_mat_is_real(A) && acb_mat_is_real(B);

    acb_mat_init(I, n, n);
    acb_mat_init(R, n, n);

    acb_mat_one(I);
    result = acb_mat_approx_solve(R, A, I, prec);
    if (result)
    {
        acb_mat_t RA, RB, E;
        mag_t d;

        acb_mat_init(RA, n, n);
        acb_mat_init(RB, n, m);
        acb_mat_init(E, n, n);
        mag_init(d);

        acb_mat_mul(RA, R, A, prec);
        acb_mat_mul(RB, R, B, prec);
        acb_mat_sub(E, I, RA, prec);
        acb_mat_bound_inf_norm(d, E);

        if (mag_cmp_2exp_si(d, 0) < 0)
        {
            slong i, j;
            mag_t e, err;

            mag_init(e);
            mag_init(err);

            mag_geom_series(d, d, 1);
            acb_mat_set(X, RB);

            for (j = 0; j < m; j++)
            {
                mag_zero(err);
                for (i = 0; i < n; i++)
                {
                    acb_get_mag(e, acb_mat_entry(RB, i, j));
                    mag_max(err, err, e);
                }
                mag_mul(err, err, d);
                for (i = 0; i < n; i++)
                {
                    if (real)
                        arb_add_error_mag(acb_realref(acb_mat_entry(X, i, j)), err);
                    else
                        acb_add_error_mag(acb_mat_entry(X, i, j), err);
                }
            }

            mag_clear(e);
            mag_clear(err);
        }
        else
        {
            result = acb_mat_solve_lu(X, RA, RB, prec);
        }

        acb_mat_clear(RA);
        acb_mat_clear(RB);
        acb_mat_clear(E);
        mag_clear(d);
    }

    acb_mat_clear(I);
    acb_mat_clear(R);

    return result;
}

int
acb_mat_solve_precond(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong n = acb_mat_nrows(A);
    slong m = acb_mat_ncols(B);

    if (m < 0.1 * n + 1)
        return _acb_mat_solve_c(X, A, B, prec);
    else
        return _acb_mat_solve_d(X, A, B, prec);
}

