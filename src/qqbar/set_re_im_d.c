/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int
qqbar_set_re_im_d(qqbar_t res, double x, double y)
{
    int ok;

    if (y == 0.0)
    {
        ok = qqbar_set_d(res, x);
    }
    else
    {
        ok = qqbar_set_d(res, y);

        if (ok)
        {
            qqbar_t t;
            qqbar_init(t);

            qqbar_i(t);
            qqbar_mul(res, res, t);

            if (x != 0)
            {
                ok = qqbar_set_d(t, x);
                qqbar_add(res, res, t);
            }

            qqbar_clear(t);
        }
    }

    return ok;
}

