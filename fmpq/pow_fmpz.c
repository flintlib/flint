/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e)
{
    slong r = WORD(1);
    FLINT_ASSERT(fmpz_sgn(e) >= 0);
    if (fmpq_is_zero(b))
    {
        if (!fmpz_is_zero(e))
            r = 0;
    }
    else if (fmpz_is_one(fmpq_denref(b)) && fmpz_is_pm1(fmpq_numref(b)))
    {
        if (!fmpz_is_one(fmpq_numref(b)) && !fmpz_is_even(e))
            r = -WORD(1);
    }
    else
    {
        if (fmpz_fits_si(e))
            fmpq_pow_si(a, b, fmpz_get_si(e));
        else
            flint_throw(FLINT_ERROR, "Exponent too large in fmpq_pow_fmpz");
        return;
    }
    fmpq_set_si(a, r, UWORD(1));
}
