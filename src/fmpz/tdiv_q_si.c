/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
{
    fmpz c1 = *g;
    slong c2 = h;

    if (h == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_tdiv_q_si). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        fmpz_set_si(f, c1 / c2);
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);

        if (c2 > 0)
        {
            flint_mpz_tdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
        }
        else
        {
            flint_mpz_tdiv_q_ui(mf, COEFF_TO_PTR(c1), -(ulong) c2);
            mpz_neg(mf, mf);
        }
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}
