/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void fmpz_get_signed_uiui(mp_limb_t * hi, mp_limb_t * lo, const fmpz_t x)
{
    ulong r0, r1, s;

    if (!COEFF_IS_MPZ(*x))
    {
        r0 = *x;
        r1 = FLINT_SIGN_EXT(r0);
    }
    else
    {
        __mpz_struct * p = COEFF_TO_PTR(*x);
        s = -(ulong)(p->_mp_size < 0);
        r0 = p->_mp_d[0];
        if (p->_mp_size > 1 || p->_mp_size < -1)
            r1 = p->_mp_d[1];
        else
            r1 = 0;

        sub_ddmmss(r1, r0, r1^s, r0^s, s, s);
    }

    *lo = r0;
    *hi = r1;
}

