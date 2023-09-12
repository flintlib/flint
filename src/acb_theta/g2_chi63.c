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
acb_theta_g2_chi63(acb_ptr r, acb_srcptr dth, slong prec)
{
    slong g = 2;
    slong n = 1 << (2 * g);
    slong orders[2] = {1, 0};
    slong i1 = acb_theta_jet_index(orders, g); /* 0 or 1 */
    slong nb = acb_theta_jet_nb(1, g + 1);
    acb_poly_t res, aux;
    acb_t t;
    ulong ab;
    slong k;

    acb_poly_init(res);
    acb_poly_init(aux);
    acb_init(t);

    acb_poly_one(res);
    for (ab = 0; ab < n; ab++)
    {
        if (!acb_theta_char_is_even(ab, g))
        {
            acb_poly_set_coeff_acb(aux, 1, &dth[nb * ab + 1 + i1]);
            acb_poly_set_coeff_acb(aux, 0, &dth[nb * ab + 1 + (1 - i1)]);
            acb_poly_mul(res, res, aux, prec);
        }
    }
    acb_poly_scalar_mul_2exp_si(res, res, -6);
    acb_const_pi(t, prec);
    acb_pow_ui(t, t, 6, prec);
    acb_poly_scalar_div(res, res, t, prec);

    for (k = 0; k <= 6; k++)
    {
        acb_poly_get_coeff_acb(&r[k], res, k);
    }

    acb_poly_clear(res);
    acb_poly_clear(aux);
    acb_clear(t);
}
