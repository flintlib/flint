/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "arb_mat.h"

int arb_mat_spd_is_lll_reduced(const arb_mat_t A, slong tol_exp, slong prec)
{
    slong g = arb_mat_nrows(A);
    fmpz_lll_t fl;
    arb_mat_t B;
    fmpz_mat_t N, U;
    arb_t c;
    int res = 1;
    slong j, k;

    arb_mat_init(B, g, g);
    fmpz_mat_init(N, g, g);
    fmpz_mat_init(U, g, g);
    arb_init(c);

    /* Set B, check error bounds on coefficients */
    for (j = 0; (j < g) && res; j++)
    {
        for (k = 0; (k < g) && res; k++)
        {
            if (mag_cmp_2exp_si(arb_radref(arb_mat_entry(A, j, k)), tol_exp - 4) > 0)
            {
                res = 0;
            }
            arb_one(c);
            arb_mul_2exp_si(c, c, tol_exp);
            arb_add_si(c, c, 1, prec);
            arb_pow_ui(c, c, FLINT_MIN(j, k), prec);
            arb_mul(arb_mat_entry(B, j, k), c, arb_mat_entry(A, j, k), prec);
        }
    }

    res = res && arb_mat_spd_get_fmpz_mat(N, B, prec);
    if (res)
    {
        /* Default Flint LLL values, except Gram */
        fmpz_lll_context_init(fl, 0.99, 0.51, GRAM, EXACT);
        fmpz_mat_one(U);
        fmpz_lll(N, U, fl);
        res = fmpz_mat_is_one(U);
    }

    arb_mat_clear(B);
    fmpz_mat_clear(N);
    fmpz_mat_clear(U);
    arb_clear(c);
    return res;
}
