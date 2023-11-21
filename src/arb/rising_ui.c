/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_hypgeom.h"

void
arb_rising_ui(arb_t y, const arb_t x, ulong n, slong prec)
{
    arb_hypgeom_rising_ui(y, x, n, prec);
}

void
arb_rising(arb_t y, const arb_t x, const arb_t n, slong prec)
{
    arb_hypgeom_rising(y, x, n, prec);
}

