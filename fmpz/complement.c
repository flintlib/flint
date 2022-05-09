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

void fmpz_complement(fmpz_t r, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        slong res = ~(*f);
        fmpz_set_si(r, res);
    }
    else
    {
        if (r != f) /* not aliased */
        {
            mpz_mock_ptr mr, mf;
            mr = _fmpz_promote(r);
            mf = COEFF_TO_PTR(*f);
            mpz_com((mpz_ptr) mr, (mpz_ptr) mf);
            _fmpz_demote_val(r);
        }
        else /* aliased */
        {
            fmpz_t tmp;
            mpz_mock_ptr mtmp, mf;
            fmpz_init(tmp);
            mtmp = _fmpz_promote(tmp);
            mf = COEFF_TO_PTR(*f);
            mpz_com((mpz_ptr) mtmp, (mpz_ptr) mf);
            _fmpz_demote_val(tmp);
            fmpz_set(r, tmp);
            fmpz_clear(tmp);
        }
    }
}
