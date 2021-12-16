/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

void fmpz_mod_poly_deflate(fmpz_mod_poly_t result, const fmpz_mod_poly_t input,
                                     ulong deflation, const fmpz_mod_ctx_t ctx)
{
    slong res_length, i;

    if (deflation == 0)
        flint_throw(FLINT_DIVZERO, "fmpz_mod_poly_deflate");

    if (input->length <= 1 || deflation == 1)
    {
        fmpz_mod_poly_set(result, input, ctx);
        return;
    }

    res_length = (input->length - 1) / deflation + 1;
    fmpz_mod_poly_fit_length(result, res_length, ctx);
    for (i = 0; i < res_length; i++)
        fmpz_set(result->coeffs + i, input->coeffs + (i * deflation));

    _fmpz_mod_poly_set_length(result, res_length);
}

