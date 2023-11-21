/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

int
acb_siegel_is_reduced(const acb_mat_t tau, slong tol_exp, slong prec)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t mat;
    acb_mat_t c;
    arb_mat_t im;
    acb_t det;
    arb_t abs, t, u;
    slong j, k;
    int res = 1;

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(c, g, g);
    arb_mat_init(im, g, g);
    acb_init(det);
    arb_init(abs);
    arb_init(t);
    arb_init(u);

    arb_one(u);
    arb_mul_2exp_si(u, u, tol_exp);

    arb_one(t);
    arb_mul_2exp_si(t, t, -1);
    arb_add(t, t, u, prec);
    for (j = 0; (j < g) && res; j++)
    {
        for (k = 0; (k < g) && res; k++)
        {
            arb_abs(abs, acb_realref(acb_mat_entry(tau, j, k)));
            if (!arb_lt(abs, t))
            {
                res = 0;
            }
        }
    }

    if (res)
    {
        acb_mat_get_imag(im, tau);
        res = arb_mat_spd_is_lll_reduced(im, tol_exp, prec);
    }

    arb_one(t);
    arb_sub(t, t, u, prec);
    for (k = 0; k < sp2gz_nb_fundamental(g); k++)
    {
        sp2gz_fundamental(mat, k);
        acb_siegel_cocycle(c, mat, tau, prec);
        acb_mat_det(det, c, prec);
        acb_abs(abs, det, prec);
        if (!arb_gt(abs, t))
        {
            res = 0;
        }
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(c);
    arb_mat_clear(im);
    acb_clear(det);
    arb_clear(abs);
    arb_clear(t);
    arb_clear(u);
    return res;
}
