/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_g2_subst_covariant(acb_poly_t r, const acb_poly_struct* powers,
    const acb_poly_t s, slong j, slong prec)
{
    slong i;
    acb_poly_t t;
    acb_t a;

    acb_poly_init(t);
    acb_init(a);

    acb_poly_zero(r);
    for (i = 0; i <= j; i++)
    {
        acb_poly_get_coeff_acb(a, s, i);
        acb_poly_scalar_mul(t, &powers[i], a, prec);
        acb_poly_add(r, r, t, prec);
    }

    acb_poly_clear(t);
    acb_clear(a);
}
