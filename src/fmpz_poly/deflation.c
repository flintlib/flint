/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

ulong _fmpz_poly_deflation(const fmpz* a, slong len)
{
    ulong i, coeff, deflation;

    if (len <= 1)
        return len;

    coeff = 1;
    while (fmpz_is_zero(a + coeff))
        coeff++;

    deflation = n_gcd(len - 1, coeff);

    while ((deflation > 1) && (coeff + deflation < (ulong) len))
    {
        for (i = 0; i < deflation - 1; i++)
        {
            coeff++;
            if (!fmpz_is_zero(a + coeff))
                deflation = n_gcd(coeff, deflation);
        }

        if (i == deflation - 1)
            coeff++;
    }

    return deflation;
}
