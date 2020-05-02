/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_printn(const ca_qqbar_t x, slong n)
{
    acb_t t;
    slong prec;

    n = FLINT_MAX(1, n);
    prec = n * 3.333 + 10;

    acb_init(t);
    ca_qqbar_get_acb(t, x, prec);

    acb_printn(t, n, ARB_STR_NO_RADIUS);
    acb_clear(t);
}

void
ca_qqbar_printnd(const ca_qqbar_t x, slong n)
{
    ca_qqbar_printn(x, n);
    flint_printf(" (deg %wd)", ca_qqbar_degree(x));
}

