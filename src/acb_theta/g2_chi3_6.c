/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_theta.h"

void
acb_theta_g2_chi3_6(acb_poly_t res, acb_srcptr dth, slong prec)
{
    slong g = 2;
    slong n = 1 << (2 * g);
    slong orders[2] = {1, 0};
    slong i1 = acb_theta_jet_index(orders, g) - 1; /* 0 or 1 */
    slong nb = acb_theta_jet_nb(1, g);
    acb_poly_struct * aux;
    acb_poly_t s;
    acb_t den;
    ulong ab;
    slong k;

    aux = flint_malloc(6 * sizeof(acb_poly_struct));
    acb_poly_init(s);
    acb_init(den);

    for (k = 0; k < 6; k++)
    {
        acb_poly_init(&aux[k]);
    }

    k = 0;
    for (ab = 0; ab < n; ab++)
    {
        if (!acb_theta_char_is_even(ab, g))
        {
            acb_poly_set_coeff_acb(&aux[k], 1, &dth[nb * ab + 1 + i1]);
            acb_poly_set_coeff_acb(&aux[k], 0, &dth[nb * ab + 1 + (1 - i1)]);
            k++;
        }
    }
    acb_poly_mul(res, &aux[0], &aux[1], prec);
    acb_poly_mul(res, res, &aux[2], prec);
    acb_poly_mul(s, &aux[3], &aux[4], prec);
    acb_poly_mul(s, s, &aux[5], prec);
    acb_poly_mul(res, res, s, prec);
    acb_const_pi(den, prec);
    acb_pow_ui(den, den, 6, prec);
    acb_poly_scalar_div(res, res, den, prec);
    acb_poly_scalar_mul_2exp_si(res, res, -6);

    acb_poly_clear(s);
    acb_clear(den);
    for (k = 0; k < 6; k++)
    {
        acb_poly_clear(&aux[k]);
    }
    flint_free(aux);
}
