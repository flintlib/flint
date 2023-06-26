/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_siegel_reduce(acb_mat_t res, fmpz_mat_t mat, const acb_mat_t tau,
                  slong prec)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t m;
    acb_mat_t cur;
    acb_mat_t star;
    acb_t det;
    arb_t abs;
    arb_t t;

    int stop = 0;
    slong j, j0;

    fmpz_mat_init(m, 2 * g, 2 * g);
    acb_mat_init(cur, g, g);
    acb_mat_init(star, g, g);
    acb_init(det);
    arb_init(abs);
    arb_init(t);

    fmpz_mat_one(mat);
    acb_mat_set(cur, tau);

    while (!stop)
    {
        /* Reduce real and imaginary parts */
        acb_siegel_reduce_imag(m, cur, prec);
        acb_siegel_transform(cur, m, cur, prec);
        fmpz_mat_mul(mat, m, mat);

        acb_siegel_reduce_real(m, cur, prec);
        acb_siegel_transform(cur, m, cur, prec);
        fmpz_mat_mul(mat, m, mat);

        /* Loop over fundamental matrices */
        j0 = -1;
        arb_one(t);

        for (j = 0; j < fmpz_mat_nb_siegel_fund(g); j++)
        {
            fmpz_mat_siegel_fund(m, j);
            acb_siegel_cocycle(star, m, cur, prec);
            acb_mat_det(det, star, prec);
            acb_abs(abs, det, prec);

            if (arb_lt(abs, t))
            {
                j0 = j;
                arb_set(t, abs);
            }
        }

        /* Apply fundamental matrix if found */
        if (j0 != -1)
        {
            fmpz_mat_siegel_fund(m, j0);
            fmpz_mat_mul(mat, m, mat);
            acb_siegel_transform(cur, mat, tau, prec);
        }
        else
            stop = 1;
    }

    acb_siegel_transform(res, mat, tau, prec);

    fmpz_mat_clear(m);
    acb_mat_clear(cur);
    acb_mat_clear(star);
    acb_clear(det);
    arb_clear(abs);
    arb_clear(t);
}
