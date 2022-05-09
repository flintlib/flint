/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fmpz-conversions.h"
#ifdef LONGSLONG
# define flint_mpz_cdiv_ui mpz_cdiv_ui
#else
# include "gmpcompat.h"
#endif

ulong
fmpz_cdiv_ui(const fmpz_t g, ulong h)
{
    fmpz c1 = *g;
    ulong r;

    if (h == UWORD(0))
        flint_throw(FLINT_DIVZERO, "fmpz_cdiv_ui\n");

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
        return flint_mpz_cdiv_ui((mpz_ptr) COEFF_TO_PTR(c1), h);
    }
}
