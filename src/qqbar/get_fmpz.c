/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_get_fmpz(fmpz_t res, const qqbar_t x)
{
    if (qqbar_degree(x) != 1 || !fmpz_is_one(QQBAR_COEFFS(x) + 1))
    {
        flint_printf("qqbar_get_fmpz: not an integer\n");
        flint_abort();
    }

    fmpz_neg(res, QQBAR_COEFFS(x));
}
