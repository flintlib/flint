/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bernoulli.h"
#include "arb.h"

void
_bernoulli_fmpq_ui_zeta(fmpz_t num, fmpz_t den, ulong n)
{
    slong prec;
    arb_t t;

    arith_bernoulli_number_denom(den, n);

    if (n % 2)
    {
        fmpz_set_si(num, -(n == 1));
        return;
    }

    if (n < BERNOULLI_SMALL_NUMER_LIMIT)
    {
        fmpz_set_si(num, _bernoulli_numer_small[n / 2]);
        return;
    }

    arb_init(t);

    for (prec = arith_bernoulli_number_size(n) + fmpz_bits(den) + 2; ; prec += 20)
    {
        arb_bernoulli_ui_zeta(t, n, prec);
        arb_mul_fmpz(t, t, den, prec);

        if (arb_get_unique_fmpz(num, t))
            break;

        flint_printf("warning: %wd insufficient precision for Bernoulli number %wu\n", prec, n);
    }

    arb_clear(t);
}

