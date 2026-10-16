/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
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
        nn_ptr yptr;
        nn_srcptr xptr;
        slong n;

        ARF_GET_MPN_READONLY(xptr, n, x);
        ARF_GET_MPN_WRITE(yptr, n, y);
        flint_mpn_copyi(yptr, xptr, n);
    }

    /* Important. */
    ARF_XSIZE(y) = ARF_XSIZE(x);
}

void
arb_set(arb_t y, const arb_t x)
{
    if (y != x)
    {
        _arf_set_inline(arb_midref(y), arb_midref(x));
        mag_set(arb_radref(y), arb_radref(x));
    }
}

void
arb_set_d(arb_t y, double x)
{
    arf_set_d(arb_midref(y), x);
    mag_zero(arb_radref(y));
}

void
arb_set_fmpz(arb_t y, const fmpz_t x)
{
    arf_set_fmpz(arb_midref(y), x);
    mag_zero(arb_radref(y));
}

void
arb_set_si(arb_t y, slong x)
{
    arf_set_si(arb_midref(y), x);
    mag_zero(arb_radref(y));
}

void
arb_set_ui(arb_t y, ulong x)
{
    arf_set_ui(arb_midref(y), x);
    mag_zero(arb_radref(y));
}
