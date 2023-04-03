/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

/*
    if r is the returned mpz, then x = r + sm, where sm is a signed 3 limb integer
    it must be safe to add +-COEFF_MAX^2*2^FLINT_BITS to the returned sm
    t is temp space since x is not to be modified
*/
mpz_srcptr _fmpz_mpoly_get_mpz_signed_uiuiui(ulong * sm, fmpz x, mpz_ptr t)
{
    mpz_ptr p;
    slong i, abs_size;
    ulong s;

    if (!COEFF_IS_MPZ(x))
    {
        sm[0] = x;
        sm[1] = FLINT_SIGN_EXT(x);
        sm[2] = FLINT_SIGN_EXT(x);
    }
    else
    {
        p = COEFF_TO_PTR(x);

        sm[0] = 0;
        sm[1] = 0;
        sm[2] = 0;

        s = FLINT_SIGN_EXT(p->_mp_size);
        abs_size = FLINT_ABS(p->_mp_size);

        if (abs_size > 3 || (abs_size == 3 && p->_mp_d[2] >= COEFF_MAX))
            return p;

        for (i = 0; i < abs_size; i++)
            sm[i] = p->_mp_d[i];

        sub_dddmmmsss(sm[2], sm[1], sm[0], s^sm[2], s^sm[1], s^sm[0], s, s, s);
    }

    mpz_set_ui(t, 0);
    return t;
}
