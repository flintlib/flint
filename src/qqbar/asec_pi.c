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
qqbar_asec_pi(slong * p, ulong * q, const qqbar_t x)
{
    if (qqbar_is_zero(x))
    {
        *p = 0;
        *q = 1;
        return 0;
    }
    else
    {
        int res;
        qqbar_t t;
        qqbar_init(t);
        qqbar_inv(t, x);
        res = qqbar_acos_pi(p, q, t);
        qqbar_clear(t);
        return res;
    }
}
