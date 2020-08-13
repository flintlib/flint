/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_rfac_uiui(fmpz_t r, ulong x, ulong n)
{
    if (n == 0)
    {
        fmpz_one(r);
    }
    else if (n == 1)
    {
        fmpz_set_ui(r, x);
    }
    else if (x == 0)
    {
        fmpz_zero(r);
    }
    else if (x <= COEFF_MAX)
    {
        _fmpz_rfac_ui(r, (fmpz *) &x, 0, n);
    }
    else
    {
        fmpz_t tmp;
        fmpz_init_set_ui(tmp, x);
        fmpz_rfac_ui(r, tmp, n);
        fmpz_clear(tmp);
    }
}
