/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* both inputs are small */
        {
            fmpz_set_si(f, c1 - c2);
        }
        else                    /* g is small, h is large */
        {
            mpz_ptr mpz3 = _fmpz_promote(f);  /* g is saved and h is large */
            mpz_ptr mpz2 = COEFF_TO_PTR(c2);
            if (c1 < WORD(0))
            {
                flint_mpz_add_ui(mpz3, mpz2, -c1);
                mpz_neg(mpz3, mpz3);
            }
            else
                flint_mpz_ui_sub(mpz3, c1, mpz2);
            _fmpz_demote_val(f);    /* may have cancelled */
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(c2))  /* g is large, h is small */
        {
            mpz_ptr mpz3 = _fmpz_promote(f);  /* h is saved and g is large */
            mpz_ptr mpz1 = COEFF_TO_PTR(c1);
            if (c2 < WORD(0))
                flint_mpz_add_ui(mpz3, mpz1, -c2);
            else
                flint_mpz_sub_ui(mpz3, mpz1, c2);
            _fmpz_demote_val(f);    /* may have cancelled */
        }
        else                    /* g and h are large */
        {
            mpz_ptr mpz3 = _fmpz_promote(f);  /* aliasing means f is already large */
            mpz_ptr mpz1 = COEFF_TO_PTR(c1);
            mpz_ptr mpz2 = COEFF_TO_PTR(c2);
            mpz_sub(mpz3, mpz1, mpz2);
            _fmpz_demote_val(f);    /* may have cancelled */
        }
    }
}
