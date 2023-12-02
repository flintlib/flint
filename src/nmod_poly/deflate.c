/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"

void
nmod_poly_deflate(nmod_poly_t result, const nmod_poly_t input, ulong deflation)
{
    slong res_length, i;

    if (deflation == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_deflate). Division by zero.\n");
    }

    if (input->length <= 1 || deflation == 1)
    {
        nmod_poly_set(result, input);
        return;
    }

    res_length = (input->length - 1) / deflation + 1;
    nmod_poly_fit_length(result, res_length);
    for (i = 0; i < res_length; i++)
        result->coeffs[i] = input->coeffs[i*deflation];

    result->length = res_length;
}
