/*
    Copyright (C) 2012 Thomas M. DuBuisson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "fmpz.h"

flint_bitcnt_t fmpz_popcnt(const fmpz_t c)
{
    mp_limb_t d;
    fmpz c1 = *c;

    if (!COEFF_IS_MPZ(c1))
    {
        if (c1 < 0)
            return 0;

        d = c1;
        return mpn_popcount(&d, 1);
    }
    else
    {
        __mpz_struct *t = COEFF_TO_PTR(c1);

        if (mpz_sgn(t) < 0)
            return 0;
        else
            return mpz_popcount(t);
    }
}
