/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_is_int_2exp_si(const arf_t x, slong e)
{
    if (arf_is_special(x))
        return arf_is_zero(x);
    else
    {
        fmpz_t t;
        int r;
        fmpz_init(t);
        arf_bot(t, x);
        r = fmpz_cmp_si(t, e) >= 0;
        fmpz_clear(t);
        return r;
    }
}

