/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_set_re_im(qqbar_t res, const qqbar_t x, const qqbar_t y)
{
    if (qqbar_is_zero(y))
    {
        qqbar_set(res, x);
    }
    else
    {
        qqbar_t t, u;
        qqbar_init(t);
        qqbar_init(u);

        qqbar_set(t, y);
        qqbar_i(u);
        qqbar_mul(t, t, u);
        qqbar_add(res, x, t);

        qqbar_clear(t);
        qqbar_clear(u);
    }
}

