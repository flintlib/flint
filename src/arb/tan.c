/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_tan(arb_t y, const arb_t x, slong prec)
{
    arb_t u;
    arb_init(u);
    arb_sin_cos(y, u, x, prec + 4);
    arb_div(y, y, u, prec);
    arb_clear(u);
}

