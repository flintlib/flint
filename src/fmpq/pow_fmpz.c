/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

int fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e)
{
    int e_sgn;

    if (fmpq_is_zero(b))
    {
        e_sgn = fmpz_sgn(e);

        if (e_sgn < 0)
            flint_throw(FLINT_ERROR, "Division by zero in fmpq_pow_fmpz");

        fmpz_set_si(fmpq_numref(a), e_sgn == 0 ? 1 : 0);
        fmpz_one(fmpq_denref(a));
        return 1;
    }
    else if (fmpz_is_one(fmpq_denref(b)) && fmpz_is_pm1(fmpq_numref(b)))
    {
        fmpz_set_si(fmpq_numref(a),
                      fmpz_is_one(fmpq_numref(b)) || fmpz_is_even(e) ? 1 : -1);
        fmpz_one(fmpq_denref(a));
        return 1;
    }
    else
    {
        if (!fmpz_fits_si(e))
            return 0;

        fmpq_pow_si(a, b, fmpz_get_si(e));
        return 1;
    }
}
