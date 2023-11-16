/*
    Copyright (C) 2018 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
_acb_dirichlet_xi(acb_t res, const acb_t s, slong prec)
{
    acb_t pi, z1, z2, z3;

    acb_init(pi);
    acb_init(z1);
    acb_init(z2);
    acb_init(z3);

    /* (s-1) * pi^(-s/2) * gamma(1 + s/2) * zeta(s) */
    acb_const_pi(pi, prec);
    acb_mul_2exp_si(z1, s, -1);
    acb_neg(z1, z1);
    acb_pow(z1, pi, z1, prec);
    acb_mul_2exp_si(z2, s, -1);
    acb_add_ui(z2, z2, 1, prec);
    acb_gamma(z2, z2, prec);
    acb_zeta(z3, s, prec);
    acb_sub_ui(res, s, 1, prec);
    acb_mul(res, res, z1, prec);
    acb_mul(res, res, z2, prec);
    acb_mul(res, res, z3, prec);

    acb_clear(pi);
    acb_clear(z1);
    acb_clear(z2);
    acb_clear(z3);
}

void acb_dirichlet_xi(acb_t res, const acb_t s, slong prec)
{
    if (!acb_is_finite(s))
    {
        acb_indeterminate(res);
    }
    else if (acb_is_one(s))
    {
        acb_one(res);
        acb_mul_2exp_si(res, res, -1);
    }
    else if ((arf_sgn(arb_midref(acb_realref(s))) < 0 &&
        !acb_contains_zero(s)) ||
        (arb_contains_si(acb_realref(s), 1) && /* also intervals around s = 1 */
        arb_contains_zero(acb_imagref(s))))
    {
        acb_sub_ui(res, s, 1, prec);
        acb_neg(res, res);
        _acb_dirichlet_xi(res, res, prec);
    }
    else
    {
        _acb_dirichlet_xi(res, s, prec);
    }
}

