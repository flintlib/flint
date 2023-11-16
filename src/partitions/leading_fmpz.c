/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "partitions.h"

void
partitions_leading_fmpz(arb_t res, const fmpz_t n, slong prec)
{
    arb_t t;
    fmpz_t c;
    slong eprec;

    arb_init(t);
    fmpz_init(c);

    eprec = prec + fmpz_bits(n) / 2;

    /* c = 24n-1 */
    fmpz_mul_ui(c, n, 24);
    fmpz_sub_ui(c, c, 1);

    /* t = pi sqrt(24n-1) / 6 */
    arb_sqrt_fmpz(t, c, eprec);
    arb_const_pi(res, eprec);
    arb_mul(t, t, res, eprec);
    arb_div_ui(t, t, 6, eprec);

    /* res = exp(t) - exp(t) / t */
    /* todo: eprec only needed for argument reduction */
    arb_exp(res, t, eprec);
    arb_div(t, res, t, prec);
    arb_sub(res, res, t, prec);

    arb_sqrt_ui(t, 12, prec);
    arb_mul(res, res, t, prec);
    arb_div_fmpz(res, res, c, prec);

    arb_clear(t);
    fmpz_clear(c);
}

