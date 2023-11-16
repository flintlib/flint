/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

int
_arf_increment_fast(arf_t x, slong prec)
{
    if (arf_sgn(x) > 0)
    {
        mp_limb_t hi, v, cy;
        mp_ptr xptr;
        mp_size_t xn;
        slong xexp;

        xexp = ARF_EXP(x);

        if (xexp >= 1 && xexp <= FLINT_BITS - 1)
        {
            ARF_GET_MPN_READONLY(xptr, xn, x);

            hi = xptr[xn - 1];
            v = hi + (UWORD(1) << (FLINT_BITS - xexp));
            cy = v < hi;

            if (cy == 0)
            {
                xptr[xn - 1] = v;
                return 0;
            }
        }
    }

    return arf_add_ui(x, x, 1, prec, ARF_RND_DOWN);
}

void
_arb_increment_fast(arb_t x, slong prec)
{
    if (_arf_increment_fast(arb_midref(x), prec))
        arf_mag_add_ulp(arb_radref(x), arb_radref(x), arb_midref(x), prec);
}

void
arb_hypgeom_rising_ui_forward(arb_t res, const arb_t x, ulong n, slong prec)
{
    arb_t t;
    ulong k;
    slong wp;

    if (n <= 1)
    {
        if (n == 0)
            arb_one(res);
        else
            arb_set_round(res, x, prec);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(n);

    arb_init(t);

    arb_add_ui(t, x, 1, wp);
    arb_mul(res, x, t, (n == 2) ? prec : wp);

    for (k = 2; k < n; k++)
    {
        _arb_increment_fast(t, wp);
        arb_mul(res, res, t, k == (n - 1) ? prec : wp);
    }

    arb_clear(t);
}

