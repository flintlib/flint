/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

slong
arf_get_si(const arf_t x, arf_rnd_t rnd)
{
    fmpz_t t;
    slong v;
    fmpz_init(t);
    arf_get_fmpz(t, x, rnd);
    if (!fmpz_fits_si(t))
    {
        flint_throw(FLINT_ERROR, "arf_get_si: result does not fit in a signed slong\n");
    }
    v = fmpz_get_si(t);
    fmpz_clear(t);
    return v;
}

