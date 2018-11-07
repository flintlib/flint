/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void fmpz_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e)
{
    slong r = WORD(1);
    FLINT_ASSERT(fmpz_sgn(e) >= 0);

    if (fmpz_is_zero(b))
    {
        if (!fmpz_is_zero(e))
            r = 0;
    }
    else if (fmpz_is_pm1(b))
    {
        if (!fmpz_is_one(b) && !fmpz_is_even(e))
            r = -WORD(1);
    }
    else
    {
        if (fmpz_abs_fits_ui(e))
            fmpz_pow_ui(a, b, fmpz_get_ui(e));
        else
            flint_throw(FLINT_ERROR, "Exponent too large in fmpz_pow_fmpz");
        return;
    }
    fmpz_set_si(a, r);
}
