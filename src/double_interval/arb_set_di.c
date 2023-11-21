/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_interval.h"
#include "arb.h"

void
arb_set_di(arb_t res, di_t x, slong prec)
{
    arf_t t, u;
    arf_init(t);
    arf_init(u);
    arf_set_d(t, x.a);
    arf_set_d(u, x.b);
    arb_set_interval_arf(res, t, u, prec);
    arf_clear(t);
    arf_clear(u);
}
