/*
    Copyright (C) 2009 William Hart

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

void
fmpz_abs(fmpz_t f1, const fmpz_t f2)
{
    if (!COEFF_IS_MPZ(*f2))  /* coeff is small */
    {
        fmpz t = FLINT_ABS(*f2);  /* Need to save value in case of aliasing */

        _fmpz_demote(f1);

        *f1 = t;
    }
    else  /* coeff is large */
    {
        /* No need to retain value in promotion, as if aliased, both already large */
        __mpz_struct * mf1 = _fmpz_promote(f1);
        mpz_abs(mf1, COEFF_TO_PTR(*f2));
    }
}
