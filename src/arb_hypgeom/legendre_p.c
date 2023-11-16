/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
arb_hypgeom_legendre_p(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, slong prec)
{
    if (arb_is_zero(m) && arb_is_int(n) && arf_sgn(arb_midref(n)) >= 0 &&
        arf_cmpabs_2exp_si(arb_midref(n), FLINT_BITS - 1) < 0)
    {
        arb_hypgeom_legendre_p_ui(res, NULL,
            arf_get_si(arb_midref(n), ARF_RND_DOWN), z, prec);
    }
    else
    {
        acb_t t, u, v;
        acb_init(t);
        acb_init(u);
        acb_init(v);
        arb_set(acb_realref(t), n);
        arb_set(acb_realref(u), m);
        arb_set(acb_realref(v), z);
        acb_hypgeom_legendre_p(t, t, u, v, type, prec);
        if (acb_is_finite(t) && acb_is_real(t))
            arb_swap(res, acb_realref(t));
        else
            arb_indeterminate(res);
        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
    }
}


