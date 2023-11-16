/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
_acb_vec_set_powers(acb_ptr xs, const acb_t x, slong len, slong prec)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
            acb_one(xs + i);
        else if (i == 1)
            acb_set_round(xs + i, x, prec);
        else if (i % 2 == 0)
            acb_mul(xs + i, xs + i / 2, xs + i / 2, prec);
        else
            acb_mul(xs + i, xs + i - 1, x, prec);
    }
}
