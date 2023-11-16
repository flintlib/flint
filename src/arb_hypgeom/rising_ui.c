/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void
arb_hypgeom_rising_ui(arb_t y, const arb_t x, ulong n, slong prec)
{
    if (n < FLINT_MAX(prec, 100))
    {
        arb_hypgeom_rising_ui_rec(y, x, n, prec);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_add_ui(t, x, n, prec);
        arb_gamma(t, t, prec);
        arb_rgamma(y, x, prec);
        arb_mul(y, y, t, prec);
        arb_clear(t);
    }
}

void
arb_hypgeom_rising(arb_t y, const arb_t x, const arb_t n, slong prec)
{
    if (arb_is_int(n) && arf_sgn(arb_midref(n)) >= 0 &&
        arf_cmpabs_ui(arb_midref(n), FLINT_MAX(prec, 100)) < 0)
    {
        arb_hypgeom_rising_ui_rec(y, x,
            arf_get_si(arb_midref(n), ARF_RND_DOWN), prec);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_add(t, x, n, prec);
        arb_gamma(t, t, prec);
        arb_rgamma(y, x, prec);
        arb_mul(y, y, t, prec);
        arb_clear(t);
    }
}

