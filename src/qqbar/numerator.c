/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_numerator(qqbar_t res, const qqbar_t y)
{
    if (qqbar_is_algebraic_integer(y))
    {
        qqbar_set(res, y);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        qqbar_denominator(t, y);
        qqbar_mul_fmpz(res, y, t);
        fmpz_clear(t);
    }
}

