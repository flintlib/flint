/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "qqbar.h"

int
qqbar_set_d(qqbar_t res, double x)
{
    arf_t t;
    fmpq_t u;
    int ok;

    arf_init(t);
    arf_set_d(t, x);

    if (arf_is_finite(t))
    {
        fmpq_init(u);
        arf_get_fmpq(u, t);
        qqbar_set_fmpq(res, u);
        ok = 1;
        fmpq_clear(u);
    }
    else
    {
        ok = 0;
    }

    arf_clear(t);
    return ok;
}
