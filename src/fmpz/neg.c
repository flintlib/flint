/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"

void
fmpz_neg(fmpz_t f1, const fmpz_t f2)
{
    if (!COEFF_IS_MPZ(*f2))     /* coeff is small */
    {
        fmpz t = -*f2;          /* Need to save value in case of aliasing */
        _fmpz_demote(f1);
        *f1 = t;
    }
    else                        /* coeff is large */
    {
        /* No need to retain value in promotion, as if aliased, both already large */
        mpz_ptr mpz_res = _fmpz_promote(f1);
        mpz_neg(mpz_res, COEFF_TO_PTR(*f2));
    }
}
