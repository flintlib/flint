/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

int fmpz_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e)
{
    int e_sgn = fmpz_sgn(e);

    if (e_sgn < 0)
    {
        flint_throw(FLINT_ERROR, "Negative exponent in fmpz_pow_fmpz");
    }
    else if (e_sgn == 0)
    {
        fmpz_one(a);
    }
    else if (fmpz_is_zero(b))
    {
        fmpz_zero(a);
    }
    else if (fmpz_is_pm1(b))
    {
        fmpz_set_si(a, fmpz_is_one(b) || fmpz_is_even(e) ? 1 : -1);
    }
    else
    {
        if (!fmpz_fits_si(e))
            return 0;

        fmpz_pow_ui(a, b, fmpz_get_si(e));
    }
    return 1;
}
