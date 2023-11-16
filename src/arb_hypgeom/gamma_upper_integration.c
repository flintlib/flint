/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void
arb_hypgeom_gamma_upper_integration(arb_t res, const arb_t s,
        const arb_t z, int regularized, slong prec)
{
    arb_t t, u;

    arb_init(t);
    arb_init(u);

    arb_one(t);
    arb_add_ui(u, s, 1, prec);

    arb_hypgeom_u_integration(u, t, u, z, prec);

    if (arb_is_finite(u))
    {
        if (regularized != 2)
        {
            arb_pow(t, z, s, prec);
            arb_mul(u, u, t, prec);

            if (regularized == 1)
            {
                arb_rgamma(t, s, prec);
                arb_mul(u, u, t, prec);
            }
        }

        arb_neg(t, z);
        arb_exp(t, t, prec);
        arb_mul(res, t, u, prec);
    }
    else
    {
        arb_indeterminate(res);
    }

    arb_clear(t);
    arb_clear(u);
}
