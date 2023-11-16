/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_cube(acb_t r, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)
    if (arb_is_zero(b))
    {
        arb_pow_ui(acb_realref(r), a, 3, prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(a))
    {
        arb_pow_ui(acb_imagref(r), b, 3, prec);
        arb_neg(acb_imagref(r), acb_imagref(r));
        arb_zero(acb_realref(r));
    }
    else
    {
        arb_t t, u, v;

        arb_init(t);
        arb_init(u);
        arb_init(v);

        arb_mul(t, a, a, prec);
        arb_mul(u, b, b, prec);
        arb_set(v, t);

        /* t = a^2 - 3b^2 */
        arb_submul_ui(t, u, 3, prec);

        /* u = -(b^2 - 3a^2) */
        arb_submul_ui(u, v, 3, prec);
        arb_neg(u, u);

        arb_mul(acb_realref(r), t, a, prec);
        arb_mul(acb_imagref(r), u, b, prec);

        arb_clear(t);
        arb_clear(u);
        arb_clear(v);
    }
#undef a
#undef b
}

