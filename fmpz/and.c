/*
    Copyright (C) 2012 Thomas M. DuBuisson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "fmpz_mini.h"
#ifdef LONGSLONG
# define flint_mpz_init_set_si mpz_init_set_si
#else
# include "gmpcompat.h"
#endif

void fmpz_and(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1, c2;
    c1 = *g;
    c2 = *h;
    if (!COEFF_IS_MPZ(c1))
    {
        if (!COEFF_IS_MPZ(c2)) /* both inputs are small */
        {
            fmpz_set_si(f, c1 & c2);
        }
        else /* g is small, h is large */
        {
            mpz_mock_t tmp;
            mpz_mock_ptr mpz3 = _fmpz_promote(f);
            flint_mpz_init_set_si((mpz_ptr) tmp, c1);
            mpz_and((mpz_ptr) mpz3, (mpz_ptr) COEFF_TO_PTR(c2), (mpz_ptr) tmp);
            _fmpz_demote_val(f);
            mpz_clear((mpz_ptr) tmp);
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(c2)) /* g is large, h is small */
        {
            mpz_mock_t tmp;
            mpz_mock_ptr mpz3 = _fmpz_promote(f);
            flint_mpz_init_set_si((mpz_ptr) tmp, c2);
            mpz_and((mpz_ptr) mpz3, (mpz_ptr) COEFF_TO_PTR(c1), (mpz_ptr) tmp);
            _fmpz_demote_val(f);
            mpz_clear((mpz_ptr) tmp);
        }
        else /* g and h are large */
        {
            mpz_mock_ptr mpz3 = _fmpz_promote(f);
            mpz_mock_ptr mpz1 = COEFF_TO_PTR(c1);
            mpz_mock_ptr mpz2 = COEFF_TO_PTR(c2);
            mpz_and((mpz_ptr) mpz3, (mpz_ptr) mpz1, (mpz_ptr) mpz2);
            _fmpz_demote_val(f);
        }
    }
}
