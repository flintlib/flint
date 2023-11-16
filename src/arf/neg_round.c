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
arf_neg_round(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        arf_neg(y, x);
        return 0;
    }
    else
    {
        int inexact;
        slong fix;
        mp_size_t xn;
        mp_srcptr xptr;

        if (y == x)
        {
            ARF_NEG(y);
            return arf_set_round(y, y, prec, rnd);
        }
        else
        {
            ARF_GET_MPN_READONLY(xptr, xn, x);
            inexact = _arf_set_round_mpn(y, &fix, xptr, xn,
                !ARF_SGNBIT(x), prec, rnd);
            _fmpz_add_fast(ARF_EXPREF(y), ARF_EXPREF(x), fix);
            return inexact;
        }
    }
}

