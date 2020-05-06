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
ca_qqbar_acsc_pi(slong * p, ulong * q, const ca_qqbar_t x)
{
    if (ca_qqbar_is_zero(x))
    {
        *p = 0;
        *q = 1;
        return 0;
    }
    else
    {
        int res;
        ca_qqbar_t t;
        ca_qqbar_init(t);
        ca_qqbar_inv(t, x);
        res = ca_qqbar_asin_pi(p, q, t);
        ca_qqbar_clear(t);
        return res;
    }
}
