/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qqbar.h"

void
qqbar_mul_2exp_si(qqbar_t res, const qqbar_t x, slong exp)
{
    slong i, d, g, h;
    fmpz * coeffs;

    d = qqbar_degree(x);

    if (qqbar_is_zero(x) || exp == 0)
    {
        qqbar_set(res, x);
        return;
    }

    if (FLINT_BIT_COUNT(d) + FLINT_BIT_COUNT((ulong) FLINT_ABS(exp)) > FLINT_BITS - 8)
    {
        flint_throw(FLINT_ERROR, "qqbar_mul_2exp_si: ludicrously large coefficients\n");
    }

    fmpz_poly_set(QQBAR_POLY(res), QQBAR_POLY(x));
    acb_mul_2exp_si(QQBAR_ENCLOSURE(res), QQBAR_ENCLOSURE(x), exp);

    coeffs = QQBAR_COEFFS(res);

    /* todo: compute valuations in advance */
    if (exp >= 0)
    {
        for (i = 1; i <= d; i++)
            fmpz_mul_2exp(coeffs + d - i, coeffs + d - i, i * exp);
    }
    else
    {
        for (i = 1; i <= d; i++)
            fmpz_mul_2exp(coeffs + i, coeffs + i, i * (-exp));
    }

    g = fmpz_val2(coeffs);
    for (i = 1; i <= d; i++)
    {
        if (!fmpz_is_zero(coeffs + i))
        {
            h = fmpz_val2(coeffs + i);
            g = FLINT_MIN(g, h);
        }
    }

    if (g != 0)
        fmpz_poly_scalar_tdiv_2exp(QQBAR_POLY(res), QQBAR_POLY(res), g);
}

