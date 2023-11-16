/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_zeta_rs_d_coeffs(arb_ptr d, const arb_t sigma, slong k, slong prec)
{
    slong j, r, m;

    arb_t u;
    arb_init(u);

    arb_one(u);
    arb_submul_ui(u, sigma, 2, prec);

    if (k == 0)
    {
        arb_one(d);
        arb_zero(d + 1);
        return;
    }

    for (j = (3 * k) / 2; j >= 0; j--)
    {
        m = 3 * k - 2 * j;

        if (m != 0)
        {
            arb_mul_2exp_si(d + j, d + j, -1);

            if (j >= 1)
                arb_addmul(d + j, d + j - 1, u, prec);

            arb_div_ui(d + j, d + j, 2 * m, prec);

            if (j >= 2)
                arb_submul_ui(d + j, d + j - 2, m + 1, prec);
        }
    }

    if (k % 2 == 0)
    {
        j = (3 * k) / 2;
        arb_zero(d + j);
        arb_set_ui(u, 2);

        for (r = j - 1; r >= 0; r--)
        {
            if ((j - r) % 2 == 0)
                arb_submul(d + j, d + r, u, prec);
            else
                arb_addmul(d + j, d + r, u, prec);

            arb_mul_ui(u, u, 4 * j - 4 * r + 2, prec);
        }
    }

    arb_zero(d + (3 * k) / 2 + 1);

    arb_clear(u);
}

