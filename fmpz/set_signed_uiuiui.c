/*
    Copyright (C) 2017 Fredrik Johansson

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
fmpz_set_signed_uiuiui(fmpz_t r, ulong hi, ulong mid, ulong lo)
{
    int negate = 0;

    if ((slong) hi < 0)
    {
        hi = -hi - ((lo != 0) || (mid != 0));
        mid = -mid - (lo != 0);
        lo = -lo;
        negate = 1;
    }

    if (hi == 0)
    {
        if (negate)
            fmpz_neg_uiui(r, mid, lo);
        else
            fmpz_set_uiui(r, mid, lo);
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(r);
        if (z->_mp_alloc < 3)
            mpz_realloc2(z, 3 * FLINT_BITS);
        z->_mp_d[0] = lo;
        z->_mp_d[1] = mid;
        z->_mp_d[2] = hi;
        z->_mp_size = negate ? -3 : 3;
    }
}

