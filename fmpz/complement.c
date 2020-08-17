/*
    Copyright (C) 2012 Thomas M. DuBuisson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void fmpz_complement(fmpz_t r, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f)) /* f is small */
    {
	slong res = ~(*f);
	fmpz_set_si(r, res);
    } else /* f is big */
    {
        if(r != f) { /* not aliased */
            __mpz_struct *ptr, *ptr2;
            ptr = _fmpz_promote(r);
            ptr2 = COEFF_TO_PTR(*f);
            mpz_com(ptr, ptr2);
            _fmpz_demote_val(r);
        } else { /* alaised */
            fmpz_t tmp;
            __mpz_struct *ptr, *ptr2;
            fmpz_init(tmp);
            ptr = _fmpz_promote(tmp);
            ptr2 = COEFF_TO_PTR(*f);
            mpz_com(ptr, ptr2);
            _fmpz_demote_val(tmp);
            fmpz_set(r,tmp);
            fmpz_clear(tmp);
        }
    }
}

