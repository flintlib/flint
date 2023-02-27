/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"

void
fmpz_lll_randtest(fmpz_lll_t fl, flint_rand_t state)
{
    double randd, delta, eta;
    rep_type rt;
    gram_type gt;

    randd = d_randtest(state);
    if (randd > 0.5 && n_randint(state, 1))
    {
        delta = 0.25 + (randd - 0.5) * 0.75;
        if (n_randint(state, 1))
            eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
        else
            eta = 0.5 + randd * (sqrt(delta) - 0.5);
    }
    else
    {
        delta = 0.25 + randd * 0.75;
        if (n_randint(state, 1))
            eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
        else
            eta = 0.5 + randd * (sqrt(delta) - 0.5);
    }
    rt = n_randint(state, 2);
    gt = n_randint(state, 2);
    fmpz_lll_context_init(fl, delta, eta, rt, gt);
}
