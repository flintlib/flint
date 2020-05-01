/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

int
ca_qqbar_set_re_im_d(ca_qqbar_t res, double x, double y)
{
    int ok;

    if (y == 0.0)
    {
        ok = ca_qqbar_set_d(res, x);
    }
    else
    {
        ok = ca_qqbar_set_d(res, y);

        if (ok)
        {
            ca_qqbar_t t;
            ca_qqbar_init(t);

            ca_qqbar_i(t);
            ca_qqbar_mul(res, res, t);

            if (x != 0)
            {
                ok = ca_qqbar_set_d(t, x);
                ca_qqbar_add(res, res, t);
            }

            ca_qqbar_clear(t);
        }
    }

    return ok;
}

