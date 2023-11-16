/*
    Copyright (C) 2018 arbguest

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

/*
 * Helper function to compute a lower bound of 1 - inf_norm(I - A*B).
 * Returns zero when this lower bound is zero.
 */
int _mag_err_complement(mag_t m,
    const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong n;
    arb_mat_t I, E;
    mag_t err;

    n = arb_mat_nrows(A);

    arb_mat_init(I, n, n);
    arb_mat_one(I);

    arb_mat_init(E, n, n);
    arb_mat_mul(E, A, B, prec);
    arb_mat_sub(E, I, E, prec);

    mag_init(err);
    arb_mat_bound_inf_norm(err, E);

    mag_one(m);
    mag_sub_lower(m, m, err);

    mag_clear(err);
    arb_mat_clear(I);
    arb_mat_clear(E);

    return !mag_is_zero(m);
}

int arb_mat_solve_preapprox(arb_mat_t X, const arb_mat_t A,
    const arb_mat_t B, const arb_mat_t R, const arb_mat_t T, slong prec)
{
    int result;
    slong m, n;
    mag_t d;

    result = 0;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    /* Use Theorem 10.2 of Rump in Acta Numerica 2010 */
    mag_init(d);
    if (_mag_err_complement(d, R, A, prec))
    {
        arb_mat_t C;

        arb_mat_init(C, n, m);
        arb_mat_mul(C, A, T, prec);
        arb_mat_sub(C, C, B, prec);
        arb_mat_mul(C, R, C, prec);

        /* Each column gets its own error bound. */
        arb_mat_set(X, T);
        {
            slong i, j;
            mag_t e, err;

            mag_init(e);
            mag_init(err);

            for (j = 0; j < m; j++)
            {
                mag_zero(err);
                for (i = 0; i < n; i++)
                {
                    arb_get_mag(e, arb_mat_entry(C, i, j));
                    mag_max(err, err, e);
                }
                mag_div(err, err, d);
                for (i = 0; i < n; i++)
                {
                    arb_add_error_mag(arb_mat_entry(X, i, j), err);
                }
            }

            mag_clear(e);
            mag_clear(err);
        }

        arb_mat_clear(C);

        result = 1;
    }

    mag_clear(d);

    return result;
}
