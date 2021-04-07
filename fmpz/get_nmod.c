/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

mp_limb_t fmpz_get_nmod(const fmpz_t aa, nmod_t mod)
{
    fmpz A = *aa;
    mp_limb_t r, SA, UA;

    if (!COEFF_IS_MPZ(A))
    {
        SA = FLINT_SIGN_EXT(A);
        UA = FLINT_ABS(A);
        NMOD_RED(r, UA, mod);
    }
    else
    {
        mpz_srcptr a = COEFF_TO_PTR(A);
        mp_srcptr ad = a->_mp_d;
        slong an = a->_mp_size;

        if (an < 0)
        {
            SA = -UWORD(1);
            an = -an;
        }
        else
        {
            SA = 0;
        }

        if (an < 5)
        {
            r = 0;
            while (an > 0)
            {
                NMOD_RED2(r, r, ad[an - 1], mod);
                an--;
            }
        }
        else
        {
            r = mpn_mod_1(ad, an, mod.n);
        }
    }

    return (SA == 0 || r == 0) ? r : (mod.n - r);
}

