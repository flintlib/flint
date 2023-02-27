/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

int fmpz_mod_poly_is_canonical(const fmpz_mod_poly_t A, const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (A->length < 0)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mod_is_canonical(A->coeffs + i, ctx))
            return 0;

        if (fmpz_is_zero(A->coeffs + i) && i + 1 == A->length)
            return 0;
    }
    return 1;
}

