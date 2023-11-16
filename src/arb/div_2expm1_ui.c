/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_div_2expm1_ui(arb_t y, const arb_t x, ulong n, slong prec)
{
    if (n < FLINT_BITS)
    {
        arb_div_ui(y, x, (UWORD(1) << n) - 1, prec);
    }
    else if (n < 1024 + prec / 32 || n > WORD_MAX / 4)
    {
        arb_t t;
        fmpz_t e;

        arb_init(t);
        fmpz_init_set_ui(e, n);

        arb_one(t);
        arb_mul_2exp_fmpz(t, t, e);
        arb_sub_ui(t, t, 1, prec);
        arb_div(y, x, t, prec);

        arb_clear(t);
        fmpz_clear(e);
    }
    else
    {
        arb_t s, t;
        slong i, b;

        arb_init(s);
        arb_init(t);

        /* x / (2^n - 1) = sum_{k>=1} x * 2^(-k*n)*/
        arb_mul_2exp_si(s, x, -n);
        arb_set(t, s);
        b = 1;

        for (i = 2; i <= prec / n + 1; i++)
        {
            arb_mul_2exp_si(t, t, -n);
            arb_add(s, s, t, prec);
            b = i;
        }

        /* error bound: sum_{k>b} x * 2^(-k*n) <= x * 2^(-b*n - (n-1)) */
        arb_mul_2exp_si(t, x, -b*n - (n-1));
        arb_abs(t, t);
        arb_add_error(s, t);

        arb_set(y, s);

        arb_clear(s);
        arb_clear(t);
    }
}

