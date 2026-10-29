/*
    Copyright (C) 2013, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"
#include "acb_dirichlet.h"

static void
arb_const_glaisher_eval(arb_t y, slong prec)
{
    acb_struct z[2];
    acb_t s;
    slong wp;

    acb_init(z + 0);
    acb_init(z + 1);
    acb_init(s);

    wp = prec + 20;

    /* zeta'(-1); acb_dirichlet_zeta_jet uses Euler-Maclaurin summation
       directly at s = -1, or the approximate functional equation at very
       high precision */
    acb_set_si(s, -1);
    acb_dirichlet_zeta_jet(z, s, 0, 2, wp);

    arb_one(y);
    arb_div_ui(y, y, 12, wp);
    arb_sub(y, y, acb_realref(z + 1), wp);
    arb_exp(y, y, wp);

    acb_clear(z + 0);
    acb_clear(z + 1);
    acb_clear(s);
}

ARB_DEF_CACHED_CONSTANT(arb_const_glaisher, arb_const_glaisher_eval)
