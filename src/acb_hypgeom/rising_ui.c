/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_rising_ui(acb_t y, const acb_t x, ulong n, slong prec)
{
    if (n < FLINT_MAX(prec, 100))
    {
        acb_hypgeom_rising_ui_rec(y, x, n, prec);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_add_ui(t, x, n, prec);
        acb_gamma(t, t, prec);
        acb_rgamma(y, x, prec);
        acb_mul(y, y, t, prec);
        acb_clear(t);
    }
}

void
acb_hypgeom_rising(acb_t y, const acb_t x, const acb_t n, slong prec)
{
    if (acb_is_int(n) && arf_sgn(arb_midref(acb_realref(n))) >= 0 &&
        arf_cmpabs_ui(arb_midref(acb_realref(n)), FLINT_MAX(prec, 100)) < 0)
    {
        acb_hypgeom_rising_ui_rec(y, x,
            arf_get_si(arb_midref(acb_realref(n)), ARF_RND_DOWN), prec);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_add(t, x, n, prec);
        acb_gamma(t, t, prec);
        acb_rgamma(y, x, prec);
        acb_mul(y, y, t, prec);
        acb_clear(t);
    }
}

