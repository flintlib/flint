/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_nint(arf_t z, const arf_t x)
{
    if (arf_is_special(x) || arf_is_int(x))
    {
        arf_set(z, x);
    }
    else
    {
        slong exp = ARF_EXP(x);

        /* now exp cannot be too large, as we would have
           caught this in arf_is_int() */
        if (COEFF_IS_MPZ(exp) || exp <= -1)
        {
            arf_zero(z);
        }
        else if (exp == 0)   /* [0.5, 1) */
        {
            if (ARF_IS_POW2(x))
                arf_zero(z);
            else
                arf_set_si(z, ARF_SGNBIT(x) ? -1 : 1);
        }
        else if (exp == 1)
        {
            nn_srcptr xp;
            slong xn, c;
            ARF_GET_MPN_READONLY(xp, xn, x);

            c = 1 + (xp[xn - 1] >= (UWORD(3) << (FLINT_BITS - 2)));
            arf_set_si(z, ARF_SGNBIT(x) ? -c : c);
        }
        else
        {
            arf_set_round(z, x, exp, ARF_RND_NEAR);
        }
    }
}

