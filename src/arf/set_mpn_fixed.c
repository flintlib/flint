/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
_arf_set_mpn_fixed(arf_t z, mp_srcptr xp, mp_size_t xn,
        mp_size_t fixn, int negative, slong prec, arf_rnd_t rnd)
{
    slong exp, exp_shift;
    int inexact;

    exp = (slong)(xn - fixn) * FLINT_BITS;

    while (xn > 0 && xp[xn-1] == 0)
    {
        xn--;
        exp -= FLINT_BITS;
    }

    if (xn == 0)
    {
        arf_zero(z);
        return 0;
    }
    else
    {
        inexact = _arf_set_round_mpn(z, &exp_shift, xp, xn, negative, prec, rnd);
        fmpz_set_si(ARF_EXPREF(z), exp + exp_shift);
        return inexact;
    }
}

