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
acb_rising_ui(acb_t y, const acb_t x, ulong n, slong prec)
{
    acb_hypgeom_rising_ui(y, x, n, prec);
}

void
acb_rising(acb_t y, const acb_t x, const acb_t n, slong prec)
{
    acb_hypgeom_rising(y, x, n, prec);
}

