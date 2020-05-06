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
ca_qqbar_log_pi_i(slong * p, ulong * q, const ca_qqbar_t x)
{
    int ok;

    ok = ca_qqbar_is_root_of_unity(p, q, x);

    if (ok)
    {
        if (*q % 2 == 0)
            *q /= 2;
        else
            *p *= 2;
    }

    while (*p > (slong) *q)
        *p -= 2 * (slong) (*q);
}
