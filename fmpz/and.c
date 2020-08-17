/*
    Copyright (C) 2012 Thomas M. DuBuisson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_and(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
        fmpz c1,c2;
        c1 = *g;
        c2 = *h;
        if(!COEFF_IS_MPZ(c1))
        {
            if(!COEFF_IS_MPZ(c2)) /* both inputs are small */
            {
                fmpz_set_si(f, c1 & c2);
            } else /* g is small, h is large */
            {
                mpz_t tmp;
                __mpz_struct *mpz3 = _fmpz_promote(f);
                flint_mpz_init_set_si(tmp, c1);
                mpz_and(mpz3, COEFF_TO_PTR(c2), tmp);
                _fmpz_demote_val(f);
                mpz_clear(tmp);
            }
        } else
        {
            if(!COEFF_IS_MPZ(c2)) /* g is large, h is small */
            {
                mpz_t tmp;
                __mpz_struct *mpz3 = _fmpz_promote(f);
                flint_mpz_init_set_si(tmp, c2);
                mpz_and(mpz3, COEFF_TO_PTR(c1), tmp);
                _fmpz_demote_val(f);
                mpz_clear(tmp);
            } else /* g and h are large */
            {
                __mpz_struct * mpz3 = _fmpz_promote(f);
                __mpz_struct * mpz1 = COEFF_TO_PTR(c1);
                __mpz_struct * mpz2 = COEFF_TO_PTR(c2);
                mpz_and(mpz3, mpz1, mpz2);
                _fmpz_demote_val(f);
            }
        }
}
