/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "fmpz.h"

ulong
fmpz_cdiv_ui(const fmpz_t g, ulong h)
{
    fmpz c1 = *g;
    ulong r;

    if (h == UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_cdiv_ui). Division by 0.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (c1 >= WORD(1))
            r = h - 1 - ((c1 - WORD(1)) % h);
        else
            r = (-c1) % h;

        return r;
    }
    else                        /* g is large */
    {
        return flint_mpz_cdiv_ui(COEFF_TO_PTR(c1), h);
    }
}
