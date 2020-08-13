/*
    Copyright (C) 2012, 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_taylor_shift(fmpz * poly, const fmpz_t c, slong len)
{
    if (len < 64)
        _fmpz_poly_taylor_shift_horner(poly, c, len);
    else
        _fmpz_poly_taylor_shift_divconquer(poly, c, len);
}

void
fmpz_poly_taylor_shift(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c)
{
    if (f != g)
        fmpz_poly_set(g, f);

    _fmpz_poly_taylor_shift(g->coeffs, c, g->length);
}

