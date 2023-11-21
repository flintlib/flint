/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"

void
arb_hurwitz_zeta(arb_t res, const arb_t s, const arb_t z, slong prec)
{
    if (!arb_contains_si(s, 1) &&
        (arb_is_positive(z) ||
            (arb_is_int(z) && arb_is_int(s) && arb_is_nonpositive(s))))
    {
        acb_t a, b, c;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_set_arb(a, s);
        acb_set_arb(b, z);
        acb_hurwitz_zeta(c, a, b, prec);
        arb_set(res, acb_realref(c));

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }
    else
    {
        arb_indeterminate(res);
    }
}

