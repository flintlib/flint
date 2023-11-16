/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
arb_contains_fmpz(const arb_t x, const fmpz_t y)
{
    int ans;
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, y);
    ans = arb_contains_arf(x, t);
    arf_clear(t);
    return ans;
}

