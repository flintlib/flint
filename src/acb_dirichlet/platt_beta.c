/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
_arb_inv_ui(arb_t res, ulong n, slong prec)
{
    arb_set_ui(res, n);
    arb_inv(res, res, prec);
}

void
acb_dirichlet_platt_beta(arb_t res, const arb_t t, slong prec)
{
    arb_t u, v;
    arb_init(u);
    arb_init(v);
    arb_log(u, t, prec);
    arb_log(v, u, prec);
    arb_div(u, v, u, prec);
    _arb_inv_ui(res, 6, prec);
    arb_add(res, res, u, prec);
    arb_clear(u);
    arb_clear(v);
}
