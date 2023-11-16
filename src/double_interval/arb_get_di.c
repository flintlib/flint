/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_interval.h"
#include "arb.h"

di_t
arb_get_di(const arb_t x)
{
    di_t res;
    if (arf_is_nan(arb_midref(x)))
    {
        res.a = -D_INF;
        res.b = D_INF;
    }
    else
    {
        arf_t t;
        arf_init(t);
        arb_get_lbound_arf(t, x, 53);
        res.a = arf_get_d(t, ARF_RND_FLOOR);
        arb_get_ubound_arf(t, x, 53);
        res.b = arf_get_d(t, ARF_RND_CEIL);
        arf_clear(t);
    }
    return res;
}
