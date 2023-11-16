/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

static inline void
_arf_set_inline(arf_t y, const arf_t x)
{
    /* Fast path */
    if (!COEFF_IS_MPZ(ARF_EXP(x)) && !COEFF_IS_MPZ(ARF_EXP(y)))
        ARF_EXP(y) = ARF_EXP(x);
    else
        fmpz_set(ARF_EXPREF(y), ARF_EXPREF(x));

    /* Fast path */
    if (!ARF_HAS_PTR(x))
    {
        ARF_DEMOTE(y);
        (y)->d = (x)->d;
    }
    else
    {
        mp_ptr yptr;
        mp_srcptr xptr;
        mp_size_t n;

        ARF_GET_MPN_READONLY(xptr, n, x);
        ARF_GET_MPN_WRITE(yptr, n, y);
        flint_mpn_copyi(yptr, xptr, n);
    }

    /* Important. */
    ARF_XSIZE(y) = ARF_XSIZE(x);
}

void
arb_set(arb_t x, const arb_t y)
{
    if (x != y)
    {
        _arf_set_inline(arb_midref(x), arb_midref(y));
        mag_set(arb_radref(x), arb_radref(y));
    }
}

