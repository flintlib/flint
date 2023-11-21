/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_floor(arf_t z, const arf_t x)
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
        if (COEFF_IS_MPZ(exp) || exp <= 0)
        {
            if (ARF_SGNBIT(x))
                arf_set_si(z, -1);
            else
                arf_zero(z);
        }
        else if (exp == 1)
        {
            arf_set_si(z, ARF_SGNBIT(x) ? -2 : 1);
        }
        else
        {
            arf_set_round(z, x, exp, ARF_RND_FLOOR);
        }
    }
}

