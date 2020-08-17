/*
    Copyright (C) 2010 Sebastian Pancratz

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
fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
{
    fmpz c1 = *g;
    slong c2 = h;

    if (h == 0)
    {
        flint_printf("Exception (fmpq_fdiv_q_si). Division by zero.\n");
        flint_abort();
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        fmpz q = c1 / c2;       /* compute C quotient */
        fmpz r = c1 - c2 * q;   /* compute remainder */

        if (r && ((c1 ^ c2) < 0))
            --q;

        fmpz_set_si(f, q);
    }
    else                        /* g is large */
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);

        if (c2 > 0)
        {
            flint_mpz_fdiv_q_ui(mpz_ptr, COEFF_TO_PTR(c1), c2);
        }
        else
        {
            flint_mpz_cdiv_q_ui(mpz_ptr, COEFF_TO_PTR(c1), -(ulong) c2);
            mpz_neg(mpz_ptr, mpz_ptr);
        }
        _fmpz_demote_val(f);    /* division by h may result in small value */
    }
}
