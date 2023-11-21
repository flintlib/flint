/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_debug(const arf_t x)
{
    mp_srcptr d;
    mp_size_t n;
    slong i;

    flint_printf("{exp=");
    fmpz_print(&x->exp);
    flint_printf("; size=%wu; sgnbit=%wd; digits=[", ARF_SIZE(x), ARF_SGNBIT(x));

    ARF_GET_MPN_READONLY(d, n, x);

    for (i = 0; i < n; i++)
        flint_printf(" %wu", d[i]);

    flint_printf("]}");
}

