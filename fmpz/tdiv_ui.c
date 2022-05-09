/*
    Copyright (C) 2009 William Hart

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
# define flint_mpz_tdiv_ui mpz_tdiv_ui
#else
# include "gmpcompat.h"
#endif

ulong
fmpz_tdiv_ui(const fmpz_t g, ulong h)
{
    fmpz c1 = *g;

    if (h == UWORD(0))
        flint_throw(FLINT_DIVZERO, "fmpz_tdiv_ui\n");

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
		/* We need the absolut value of the remainder and
		   C 90 guarantees truncation towards zero. */
		if (c1 < WORD(0))
			return -c1 % h;
		else
			return c1 % h;
    }
    else                        /* g is large */
    {
        return flint_mpz_tdiv_ui((mpz_ptr) COEFF_TO_PTR(c1), h);
    }
}
