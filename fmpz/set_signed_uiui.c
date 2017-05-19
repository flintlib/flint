/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_set_signed_uiui(fmpz_t r, ulong hi, ulong lo)
{
    if (((slong) hi) < 0)
    {
        hi = -hi - (lo != 0);
        lo = -lo;
        fmpz_neg_uiui(r, hi, lo);
    }
    else
    {
        fmpz_set_uiui(r, hi, lo);
    }
}

