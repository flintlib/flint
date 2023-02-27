/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
fmpz_poly_inflate(fmpz_poly_t result, const fmpz_poly_t input, ulong inflation)
{
    if (input->length <= 1 || inflation == 1)
    {
        fmpz_poly_set(result, input);
    }
    else if (inflation == 0)
    {
        fmpz_t v;
        fmpz_init_set_ui(v, 1);
        fmpz_poly_evaluate_fmpz(v, input, v);
        fmpz_poly_zero(result);
        fmpz_poly_set_coeff_fmpz(result, 0, v);
        fmpz_clear(v);
    }
    else
    {
        slong i, j, res_length = (input->length - 1) * inflation + 1;

        fmpz_poly_fit_length(result, res_length);

        for (i = input->length - 1; i > 0; i--)
        {
            fmpz_set(result->coeffs + (i * inflation), input->coeffs + i);
            for (j = i * inflation - 1; j > (i - 1) * inflation; j--)
                fmpz_zero(result->coeffs + j);
        }
        fmpz_set(result->coeffs + 0, input->coeffs + 0);
        result->length = res_length;
    }
}

