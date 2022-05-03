/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#ifdef LONGSLONG
# define flint_mpz_divexact_ui mpz_divexact_ui
#else
# include "gmpcompat.h"
#endif

void fmpz_divexact_si(fmpz_t f, const fmpz_t g, slong h)
{
    fmpz c1 = *g;

    if (h == 0)
        flint_throw(FLINT_DIVZERO, "fmpz_divexact_si\n");

    if (!COEFF_IS_MPZ(c1))  /* g is small */
    {
        fmpz_set_si(f, c1 / h);
    }
    else  /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);

        if (h > 0)
        {
            flint_mpz_divexact_ui(mf, COEFF_TO_PTR(c1), h);
            _fmpz_demote_val(f);  /* division by h may result in small value */
        }
        else
        {
            flint_mpz_divexact_ui(mf, COEFF_TO_PTR(c1), -((ulong) h));
            _fmpz_demote_val(f);  /* division by h may result in small value */
            fmpz_neg(f, f);
        }
    }
}
