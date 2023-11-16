/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "arb_fmpz_poly.h"

ulong arb_fmpz_poly_deflation(const fmpz_poly_t input)
{
    slong i, coeff, deflation;

    if (input->length <= 1)
        return input->length;

    coeff = 1;
    while (fmpz_is_zero(input->coeffs + coeff))
        coeff++;

    deflation = n_gcd(input->length - 1, coeff);

    while ((deflation > 1) && (coeff + deflation < input->length))
    {
        for (i = 0; i < deflation - 1; i++)
        {
            coeff++;
            if (!fmpz_is_zero(input->coeffs + coeff))
                deflation = n_gcd(coeff, deflation);
        }

        if (i == deflation - 1)
            coeff++;
    }

    return deflation;
}

