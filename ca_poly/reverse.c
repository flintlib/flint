/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_reverse(ca_ptr res, ca_srcptr poly, slong len, slong n, ca_ctx_t ctx)
{
    if (res == poly)
    {
        slong i;

        for (i = 0; i < n / 2; i++)
        {
            ca_struct t = res[i];
            res[i] = res[n - 1 - i];
            res[n - 1 - i] = t;
        }

        for (i = 0; i < n - len; i++)
            ca_zero(res + i, ctx);
    }
    else
    {
        slong i;

        for (i = 0; i < n - len; i++)
            ca_zero(res + i, ctx);

        for (i = 0; i < len; i++)
            ca_set(res + (n - len) + i, poly + (len - 1) - i, ctx);
    }
}
