/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "fmpz.h"

flint_bitcnt_t fmpz_val2(const fmpz_t x)
{
    fmpz c = *x;
    flint_bitcnt_t t;

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 0)
            t = 0;
        else
            t = flint_ctz(FLINT_ABS(c));
    }
    else
    {
        ulong *d = (COEFF_TO_PTR(c))->_mp_d;
        flint_bitcnt_t u;

        t = 0;
        while (*d == 0)
        {
            d++;
            t += FLINT_BITS;
        }

        u = flint_ctz(*d);
        t += u;
    }

    return t;
}
