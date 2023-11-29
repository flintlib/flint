/*
    Copyright (C) 2009 William Hart

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
fmpz_tdiv_ui(const fmpz_t g, ulong h)
{
    fmpz c1 = *g;

    if (h == UWORD(0))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_tdiv_ui). Division by 0.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
		/* We need the absolute value of the remainder and
		   C 90 guarantees truncation towards zero. */
		if (c1 < WORD(0))
			return -c1 % h;
		else
			return c1 % h;
    }
    else                        /* g is large */
    {
        return flint_mpz_tdiv_ui(COEFF_TO_PTR(c1), h);
    }
}
