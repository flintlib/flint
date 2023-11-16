/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_quadratic_roots_fmpz(acb_t r1, acb_t r2,
    const fmpz_t a, const fmpz_t b, const fmpz_t c, slong prec)
{
    fmpz_t d;
    fmpz_init(d);

    /* d = b^2 - 4ac */
    fmpz_mul(d, a, c);
    fmpz_mul_2exp(d, d, 2);
    fmpz_submul(d, b, b);
    fmpz_neg(d, d);

    /* +/- sqrt(d) */
    acb_zero(r1);
    if (fmpz_sgn(d) >= 0)
    {
        arb_sqrt_fmpz(acb_realref(r1), d, prec + fmpz_bits(d) + 4);
    }
    else
    {
        fmpz_neg(d, d);
        arb_sqrt_fmpz(acb_imagref(r1), d, prec + fmpz_bits(d) + 4);
    }
    acb_neg(r2, r1);

    /* -b */
    acb_sub_fmpz(r1, r1, b, prec + 4);
    acb_sub_fmpz(r2, r2, b, prec + 4);

    /* divide by 2a */
    fmpz_mul_2exp(d, a, 1);
    acb_div_fmpz(r1, r1, d, prec);
    acb_div_fmpz(r2, r2, d, prec);

    fmpz_clear(d);
    return;
}

