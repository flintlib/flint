/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "double_extras.h"
#include "fmpz_lll.h"

#define EPS 0.000000000003

void
fmpz_lll_randtest(fmpz_lll_t fl, flint_rand_t state)
{
    double delta, eta;
    rep_type rt;
    gram_type gt;

    delta = 0.25 + 0.75 * d_randtest(state);
    if (delta <= 0.25)
        delta = nextafter(0.25, 1.00);

    eta = 0.5 + 2 * (d_randtest(state) - 0.5) * (sqrt(delta) - 0.5);
    if (eta <= 0.5 + EPS)
        eta = 0.5 + 8 * EPS * d_randtest(state);

    rt = n_randint(state, 2);
    gt = n_randint(state, 2);

    fmpz_lll_context_init(fl, delta, eta, rt, gt);
}
