/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_poly.h"

void
arb_const_glaisher_eval(arb_t y, slong prec)
{
    acb_struct z[2];
    acb_t s, a;
    slong wp;

    acb_init(z + 0);
    acb_init(z + 1);
    acb_init(s);
    acb_init(a);

    wp = prec + 20;

    /* directly evaluating at s = -1 is slightly faster
       than evaluating at s = 2 */
    acb_set_si(s, -1);
    acb_one(a);
    _acb_poly_zeta_cpx_series(z, s, a, 0, 2, wp);

    arb_one(y);
    arb_div_ui(y, y, 12, wp);
    arb_sub(y, y, acb_realref(z + 1), wp);
    arb_exp(y, y, wp);

    acb_clear(z + 0);
    acb_clear(z + 1);
    acb_clear(s);
    acb_clear(a);
}

ARB_DEF_CACHED_CONSTANT(arb_const_glaisher, arb_const_glaisher_eval)

