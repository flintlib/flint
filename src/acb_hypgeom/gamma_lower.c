/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
{
    acb_t s1, nz, t, u;

    if (!acb_is_finite(s) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(s1);
    acb_init(nz);
    acb_init(t);
    acb_init(u);

    acb_add_ui(s1, s, 1, prec);
    acb_neg(nz, z);

    if (regularized == 0)
    {
        /* \gamma(s, z) = s^-1 z^s 1F1(s, 1+s, -z) */
        acb_hypgeom_m(u, s, s1, nz, 0, prec);
        acb_pow(t, z, s, prec);
        acb_mul(u, u, t, prec);
        acb_div(res, u, s, prec);
    }
    else if (regularized == 1)
    {
        /* P(s, z) = z^s \gamma^{*}(s, z) */
        acb_hypgeom_m(u, s, s1, nz, 1, prec);
        acb_pow(t, z, s, prec);
        acb_mul(res, u, t, prec);
    }
    else if (regularized == 2)
    {
        /* \gamma^{*}(s, z) */
        acb_hypgeom_m(res, s, s1, nz, 1, prec);
    }

    acb_clear(s1);
    acb_clear(nz);
    acb_clear(t);
    acb_clear(u);
}

