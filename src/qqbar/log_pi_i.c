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
qqbar_log_pi_i(slong * p, ulong * q, const qqbar_t x)
{
    int ok;

    ok = qqbar_is_root_of_unity(p, q, x);

    if (ok)
    {
        if (*q % 2 == 0)
            *q /= 2;
        else
            *p *= 2;

        while (*p > (slong) *q)
            *p -= 2 * (slong) (*q);
    }

    return ok;
}
