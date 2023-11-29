/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "fmpz.h"

void
fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
{
    fmpz c1 = *g;
    ulong c2 = h;

    if (h == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_fdiv_q_ui). Division by zero.\n");
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (c1 > 0)
        {
            fmpz_set_ui(f, c1 / c2);
        }
        else
        {
            ulong q = ((ulong) -c1) / c2;
            ulong r = ((ulong) -c1) - c2 * q;

            if (r)
                ++q;

            fmpz_set_si(f, - (slong) q);
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);

        flint_mpz_fdiv_q_ui(mf, COEFF_TO_PTR(c1), c2);
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}
