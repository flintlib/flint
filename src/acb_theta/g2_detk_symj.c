/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_g2_detk_symj(acb_poly_t res, const acb_mat_t m, const acb_poly_t f,
    slong k, slong j, slong prec)
{
    acb_poly_t x, y, t, u, aux;
    acb_t a;
    slong i;

    acb_poly_init(x);
    acb_poly_init(y);
    acb_poly_init(t);
    acb_poly_init(u);
    acb_poly_init(aux);
    acb_init(a);

    acb_poly_set_coeff_acb(x, 0, acb_mat_entry(m, 1, 0));
    acb_poly_set_coeff_acb(x, 1, acb_mat_entry(m, 0, 0));
    acb_poly_set_coeff_acb(y, 0, acb_mat_entry(m, 1, 1));
    acb_poly_set_coeff_acb(y, 1, acb_mat_entry(m, 0, 1));

    for (i = 0; i <= j; i++)
    {
        acb_poly_get_coeff_acb(a, f, i);
        acb_poly_pow_ui(t, x, i, prec);
        acb_poly_pow_ui(u, y, j - i, prec);
        acb_poly_mul(t, t, u, prec);
        acb_poly_scalar_mul(t, t, a, prec);
        acb_poly_add(aux, aux, t, prec);
    }
    acb_mat_det(a, m, prec);
    acb_pow_si(a, a, k, prec);
    acb_poly_scalar_mul(res, aux, a, prec);

    acb_poly_clear(x);
    acb_poly_clear(y);
    acb_poly_clear(aux);
    acb_poly_clear(t);
    acb_poly_clear(u);
    acb_clear(a);
}
