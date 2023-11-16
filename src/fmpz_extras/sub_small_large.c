/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_extras.h"

slong
_fmpz_sub_small_large(const fmpz_t x, const fmpz_t y)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_sub(t, x, y);
    if (!COEFF_IS_MPZ(*t))
    {
        /* no need to free t */
        return *t;
    }
    else
    {
        int sign = fmpz_sgn(t);
        fmpz_clear(t);
        return (sign > 0) ? WORD_MAX : -WORD_MAX;
    }
}
