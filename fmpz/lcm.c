/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void
fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz_t t;

    if (fmpz_is_zero(g) || fmpz_is_zero(h))
    {
        fmpz_zero(f);
        return;
    }

    if (fmpz_is_pm1(g))
    {
        fmpz_abs(f, h);
        return;
    }

    if (fmpz_is_pm1(h))
    {
        fmpz_abs(f, g);
        return;
    }

    fmpz_init(t);
    fmpz_gcd(t, g, h);
    fmpz_divexact(t, g, t);
    fmpz_mul(f, t, h);
    fmpz_abs(f, f);
    fmpz_clear(t);
}
