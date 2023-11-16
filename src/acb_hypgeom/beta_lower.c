/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

/* todo: move this to the acb module? */
static void
acb_beta(acb_t res, const acb_t a, const acb_t b, slong prec)
{
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    acb_gamma(t, a, prec);
    acb_gamma(u, b, prec);

    acb_add(res, a, b, prec);
    acb_rgamma(res, res, prec);
    acb_mul(res, res, t, prec);
    acb_mul(res, res, u, prec);

    acb_clear(t);
    acb_clear(u);
}

void acb_hypgeom_beta_lower(acb_t res,
    const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    acb_t t, u;

    if (acb_is_zero(z) && arb_is_positive(acb_realref(a)))
    {
        acb_zero(res);
        return;
    }

    if (acb_is_one(z) && arb_is_positive(acb_realref(b)))
    {
        if (regularized)
            acb_one(res);
        else
            acb_beta(res, a, b, prec);
        return;
    }

    acb_init(t);
    acb_init(u);

    acb_sub_ui(t, b, 1, prec);
    acb_neg(t, t);
    acb_add_ui(u, a, 1, prec);

    if (regularized)
    {
        acb_hypgeom_2f1(t, a, t, u, z, 1, prec);

        acb_add(u, a, b, prec);
        acb_gamma(u, u, prec);
        acb_mul(t, t, u, prec);
        acb_rgamma(u, b, prec);
        acb_mul(t, t, u, prec);
    }
    else
    {
        acb_hypgeom_2f1(t, a, t, u, z, 0, prec);
        acb_div(t, t, a, prec);
    }

    acb_pow(u, z, a, prec);
    acb_mul(t, t, u, prec);

    acb_set(res, t);

    acb_clear(t);
    acb_clear(u);
}

