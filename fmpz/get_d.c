/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz-impl.h"

double
fmpz_get_d(const fmpz_t f)
{
    fmpz c = *f;

    if (c >= DOUBLE_MIN && c <= DOUBLE_MAX)
    {
        return (double) c;
    }
    else if (!COEFF_IS_MPZ(c))
    {
        mp_limb_t d;

        if (c > 0)
        {
            d = c;
            return flint_mpn_get_d(&d, 1, 1, 0);
        }
        else
        {
            d = -c;
            return flint_mpn_get_d(&d, 1, -1, 0);
        }
    }
    else
        return mpz_get_d(COEFF_TO_PTR(c));
}
