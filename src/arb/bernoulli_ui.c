/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb.h"
#include "bernoulli.h"

void
arb_bernoulli_ui(arb_t b, ulong n, slong prec)
{
    if (n < bernoulli_cache_num)
    {
        arb_set_fmpq(b, bernoulli_cache + n, prec);
    }
    else
    {
        int use_frac;

        use_frac = (n < BERNOULLI_SMALL_NUMER_LIMIT) || (n % 2 != 0);
        if (!use_frac && n < UWORD_MAX / 1000)
            use_frac = (prec > bernoulli_global_prec(n));

        if (use_frac)
        {
            fmpq_t t;
            fmpq_init(t);
            bernoulli_fmpq_ui(t, n);
            arb_set_fmpq(b, t, prec);
            fmpq_clear(t);
        }
        else
        {
            arb_bernoulli_ui_zeta(b, n, prec);
        }
    }
}

