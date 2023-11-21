/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_reverse(arb_ptr res, arb_srcptr poly, slong len, slong n)
{
    if (res == poly)
    {
        slong i;

        for (i = 0; i < n / 2; i++)
        {
            arb_struct t = res[i];
            res[i] = res[n - 1 - i];
            res[n - 1 - i] = t;
        }

        for (i = 0; i < n - len; i++)
            arb_zero(res + i);
    }
    else
    {
        slong i;

        for (i = 0; i < n - len; i++)
            arb_zero(res + i);

        for (i = 0; i < len; i++)
            arb_set(res + (n - len) + i, poly + (len - 1) - i);
    }
}
