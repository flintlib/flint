/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

ulong
nmod_poly_deflation(const nmod_poly_t input)
{
    ulong deflation;
    slong i, coeff;

    if (input->length <= 1)
        return input->length;

    coeff = 1;
    while (!input->coeffs[coeff])
        coeff++;

    deflation = n_gcd(input->length - 1, coeff);

    while ((deflation > 1) && (coeff + deflation < input->length))
    {
        for (i = 0; i < deflation - 1; i++)
        {
            coeff++;
            if (input->coeffs[coeff])
                deflation = n_gcd(coeff, deflation);
        }
        if (i == deflation - 1)
            coeff++;
    }

    return deflation;
}
