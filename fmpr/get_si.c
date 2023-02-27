/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
fmpr_get_si(const fmpr_t x, fmpr_rnd_t rnd)
{
    fmpz_t t;
    slong v;
    fmpz_init(t);
    fmpr_get_fmpz(t, x, rnd);
    if (!fmpz_fits_si(t))
    {
        flint_printf("fmpr_get_si: result does not fit in a signed slong\n");
        flint_abort();
    }
    v = fmpz_get_si(t);
    fmpz_clear(t);
    return v;
}

