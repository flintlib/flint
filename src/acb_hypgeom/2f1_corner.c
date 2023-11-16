/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_2f1_corner(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, slong prec)
{
    acb_t aa, bb, cc, z1, z2, f1, f2;
    int upper;

    acb_init(aa); acb_init(bb); acb_init(cc);
    acb_init(z1); acb_init(z2); acb_init(f1); acb_init(f2);

    acb_add_ui(aa, a, 1, prec);
    acb_add_ui(bb, b, 1, prec);
    acb_add_ui(cc, c, 1, prec);

    upper = arb_is_positive(acb_imagref(z));

    /* 0 -> 0.5 +/- 0.5i -> 0.5 +/- 0.75i -> z */
#if 0
    acb_set_d_d(z1, 0.5, upper ? 0.5 : -0.5);
    acb_set_d_d(z2, 0.5, upper ? 0.75 : -0.75);
#else
    acb_set_d_d(z1, 0.375, upper ? 0.625 : -0.625);
    acb_set_d_d(z2, 0.5, upper ? 0.8125 : -0.8125);
#endif

    acb_hypgeom_2f1_direct(f1, a, b, c, z1, regularized, prec);

    acb_hypgeom_2f1_direct(f2, aa, bb, cc, z1, regularized, prec);
    acb_mul(f2, f2, a, prec);
    acb_mul(f2, f2, b, prec);
    if (!regularized)
        acb_div(f2, f2, c, prec);

    acb_hypgeom_2f1_continuation(f1, f2, a, b, c, z1, z2, f1, f2, prec);

    acb_set(z1, z2);
    acb_set(z2, z);

    acb_hypgeom_2f1_continuation(f1, f2, a, b, c, z1, z2, f1, f2, prec);

    acb_set(res, f1);

    acb_clear(aa); acb_clear(bb); acb_clear(cc);
    acb_clear(z1); acb_clear(z2); acb_clear(f1); acb_clear(f2);
}

