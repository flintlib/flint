/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_poly.h"

void
acb_polygamma(acb_t res, const acb_t s, const acb_t z, slong prec)
{
    if (acb_is_zero(s))
    {
        acb_digamma(res, z, prec);
    }
    else if (acb_is_int(s) && arb_is_positive(acb_realref(s)))
    {
        acb_t t, u;

        acb_init(t);
        acb_init(u);

        acb_add_ui(t, s, 1, prec);
        acb_gamma(u, t, prec);
        acb_hurwitz_zeta(t, t, z, prec);

        if (acb_is_int_2exp_si(s, 1))
            acb_neg(t, t);

        acb_mul(res, t, u, prec);

        acb_clear(t);
        acb_clear(u);
    }
    else
    {
        acb_t t, u;
        acb_struct v[2];

        acb_init(t);
        acb_init(u);

        acb_init(v);
        acb_init(v + 1);

        /* u = psi(-s) + gamma */
        acb_neg(t, s);
        acb_digamma(u, t, prec);
        arb_const_euler(acb_realref(v), prec);
        arb_add(acb_realref(u), acb_realref(u), acb_realref(v), prec);

        acb_add_ui(t, s, 1, prec);
        _acb_poly_zeta_cpx_series(v, t, z, 0, 2, prec);

        acb_addmul(v + 1, v, u, prec);

        acb_neg(t, s);
        acb_rgamma(u, t, prec);
        acb_mul(res, v + 1, u, prec);

        acb_clear(v);
        acb_clear(v + 1);

        acb_clear(t);
        acb_clear(u);
    }
}

