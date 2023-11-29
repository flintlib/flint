/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

void
fmpz_poly_deflate(fmpz_poly_t result, const fmpz_poly_t input, ulong deflation)
{
    slong res_length, i;

    if (deflation == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_deflate). Division by zero.\n");
    }

    if (input->length <= 1 || deflation == 1)
    {
        fmpz_poly_set(result, input);
        return;
    }

    res_length = (input->length - 1) / deflation + 1;
    fmpz_poly_fit_length(result, res_length);
    for (i = 0; i < res_length; i++)
        fmpz_set(result->coeffs + i, input->coeffs + i*deflation);

    result->length = res_length;
}

