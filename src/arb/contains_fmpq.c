/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
arb_contains_fmpq(const arb_t x, const fmpq_t y)
{
    if (fmpz_is_one(fmpq_denref(y)) || !arb_is_finite(x))
    {
        return arb_contains_fmpz(x, fmpq_numref(y));
    }
    else
    {
        arf_t t, xm, xr, ym;
        arf_struct tmp[3];
        int result;

        arf_init(t);
        arf_init(xm);
        arf_init(xr);
        arf_init(ym);

        /* To compare x with p/q, compare qx with p. */
        arf_mul_fmpz(xm, arb_midref(x), fmpq_denref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_mag(xr, arb_radref(x));
        arf_mul_fmpz(xr, xr, fmpq_denref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_fmpz(ym, fmpq_numref(y));

        /* y >= xm - xr  <=>  0 >= xm - xr - y */
        arf_init_set_shallow(tmp + 0, xm);
        arf_init_neg_shallow(tmp + 1, xr);
        arf_init_neg_shallow(tmp + 2, ym);

        arf_sum(t, tmp, 3, 30, ARF_RND_DOWN);
        result = (arf_sgn(t) <= 0);

        if (result)
        {
            /* y <= xm + xr  <=>  0 <= xm + xr - y */
            arf_init_set_shallow(tmp + 1, xr);
            arf_sum(t, tmp, 3, 30, ARF_RND_DOWN);
            result = (arf_sgn(t) >= 0);
        }

        arf_clear(t);
        arf_clear(xm);
        arf_clear(xr);
        arf_clear(ym);

        return result;
    }
}

