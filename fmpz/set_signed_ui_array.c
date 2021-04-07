/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

/*
    Given an array of limbs "c" representing a integer mod 2^(FLINT_BITS*n),
    set "f" to the symmetric remainder with the halfway point
    2^(FLINT_BITS*n/2) mapping to -2^(FLINT_BITS*n/2)
*/

void fmpz_set_signed_ui_array(fmpz_t f, const ulong * c, slong n)
{
    ulong csign;

    FLINT_ASSERT(n > 0);

    csign = FLINT_SIGN_EXT(c[n - 1]);

    while (n > 0 && c[n - 1] == csign)
        n--;

    if (n < 2)
    {
        if (csign == 0)
            fmpz_set_ui(f, c[0]);
        else if (c[0] != 0)
            fmpz_neg_ui(f, -c[0]);
        else
            fmpz_neg_uiui(f, 1, 0);
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(f);
        mp_limb_t * zd = FLINT_MPZ_REALLOC(z, n);

        if (csign == 0)
        {
            flint_mpn_copyi(zd, c, n);
            z->_mp_size = n;
        }
        else
        {
            if (mpn_neg(zd, c, n))
            {
                FLINT_ASSERT(zd[n - 1] != 0);
                z->_mp_size = -n;
            }
            else
            {
                zd = FLINT_MPZ_REALLOC(z, n + 1);
                zd[n] = 1;
                z->_mp_size = -(n + 1);
            }
        }
    }
}

