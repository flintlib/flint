/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"

void
fmpz_gcd_ui(fmpz_t res, const fmpz_t a, ulong b)
{
    if (b == 0)
        fmpz_abs(res, a);
    else if (!COEFF_IS_MPZ(*a))
    {
        if (*a != 0)
        {
            _fmpz_demote(res);
            *res = mpn_gcd_1(&b, 1, FLINT_ABS(*a));
        }
        else
            fmpz_set_ui(res, b);
    }
    else
    {
        __mpz_struct * ma = COEFF_TO_PTR(*a);
        fmpz_set_ui(res, mpn_gcd_1(ma->_mp_d, FLINT_ABS(ma->_mp_size), b));
    }
}
