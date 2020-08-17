/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_taylor_shift_horner(fmpz * poly, const fmpz_t c, slong n)
{
    slong i, j;

    if (*c == WORD(1))
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                fmpz_add(poly + j, poly + j, poly + j + 1);
    }
    else if (*c == WORD(-1))
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                fmpz_sub(poly + j, poly + j, poly + j + 1);
    }
    else if (*c != WORD(0))
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                fmpz_addmul(poly + j, poly + j + 1, c);
    }
}

void
fmpz_poly_taylor_shift_horner(fmpz_poly_t g, const fmpz_poly_t f,
    const fmpz_t c)
{
    if (f != g)
        fmpz_poly_set(g, f);

    _fmpz_poly_taylor_shift_horner(g->coeffs, c, g->length);
}
