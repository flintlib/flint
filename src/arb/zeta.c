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
arb_zeta(arb_t y, const arb_t s, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_set_arb(t, s),
    acb_zeta(t, t, prec);
    arb_set(y, acb_realref(t));
    acb_clear(t);
}

