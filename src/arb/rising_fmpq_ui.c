/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

static void
bsplit(arb_t y, const fmpz_t p, const fmpz_t q, ulong a, ulong b, slong prec)
{
    if (b - a <= 8)
    {
        fmpz_t t, u;
        ulong c;

        fmpz_init(t);
        fmpz_init(u);

        fmpz_mul_ui(t, q, a);
        fmpz_add(t, t, p);
        fmpz_set(u, t);

        for (c = a + 1; c < b; c++)
        {
            fmpz_add(u, u, q);
            fmpz_mul(t, t, u);
        }

        arb_set_round_fmpz(y, t, prec);

        fmpz_clear(t);
        fmpz_clear(u);
    }
    else
    {
        arb_t w;
        ulong m = a + (b - a) / 2;
        arb_init(w);

        bsplit(y, p, q, a, m, prec);
        bsplit(w, p, q, m, b, prec);

        arb_mul(y, y, w, prec);
        arb_clear(w);
    }
}

void
arb_rising_fmpq_ui(arb_t y, const fmpq_t x, ulong n, slong prec)
{
    if (n == 0)
    {
        arb_one(y);
    }
    else if (n == 1)
    {
        arb_set_fmpq(y, x, prec);
    }
    else
    {
        slong wp;

        wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        bsplit(y, fmpq_numref(x), fmpq_denref(x), 0, n, wp);

        if (fmpz_is_one(fmpq_denref(x)))
        {
            arb_set_round(y, y, prec);
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_set_fmpz(t, fmpq_denref(x));
            arb_pow_ui(t, t, n, wp);
            arb_div(y, y, t, prec);
            arb_clear(t);
        }
    }
}

