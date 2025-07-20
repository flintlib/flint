/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

static void
_fmpz_add_inline(fmpz_t z, const fmpz_t x, const fmpz_t y)
{
    fmpz f, g;

    f = *x;
    g = *y;

    if (!COEFF_IS_MPZ(f) && !COEFF_IS_MPZ(g))
        fmpz_set_si(z, f + g);
    else
        fmpz_add(z, x, y);
}

void
_fmpz_vec_add(fmpz * res, const fmpz * vec1, const fmpz * vec2, slong len2)
{
    slong i;
    for (i = 0; i < len2; i++)
        _fmpz_add_inline(res + i, vec1 + i, vec2 + i);
}
