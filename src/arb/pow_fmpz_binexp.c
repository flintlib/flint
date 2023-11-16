/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, slong prec)
{
    slong i, wp, bits;

    if (-WORD(2) <= *e && *e <= WORD(2))
    {
        if (*e == WORD(0))
            arb_one(y);
        else if (*e == WORD(1))
            arb_set_round(y, b, prec);
        else if (*e == -WORD(1))
            arb_inv(y, b, prec);
        else if (*e == WORD(2))
            arb_mul(y, b, b, prec);
        else
        {
            arb_inv(y, b, prec);
            arb_mul(y, y, y, prec);
        }
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);

        if (arb_is_exact(b))
        {
            arb_pow_fmpz_binexp(y, b, f, prec + 2);
            arb_inv(y, y, prec);
        }
        else
        {
            arb_inv(y, b, prec + fmpz_bits(e) + 2);
            arb_pow_fmpz_binexp(y, y, f, prec);
        }

        fmpz_clear(f);
        return;
    }

    if (y == b)
    {
        arb_t t;
        arb_init(t);
        arb_set(t, b);
        arb_pow_fmpz_binexp(y, t, e, prec);
        arb_clear(t);
        return;
    }

    arb_set(y, b);

    bits = fmpz_bits(e);
    wp = ARF_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        arb_mul(y, y, y, wp);
        if (fmpz_tstbit(e, i))
            arb_mul(y, y, b, wp);
    }
}

