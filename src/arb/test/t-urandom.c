/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

#define N 10000

TEST_FUNCTION_START(arb_urandom, state)
{
    slong iter;
    slong prec;
    arb_ptr rand;
    arb_t m; /* mean */
    arb_t s; /* variance */
    arb_t mp;
    arb_t sp;
    arb_t tmp;

    arb_init(m);
    arb_init(s);
    arb_init(mp);
    arb_init(sp);
    arb_init(tmp);
    rand = _arb_vec_init(N);
    prec = 299;

    for (iter = 0; iter < N; iter++)
    {
        arb_urandom(rand + iter, state, prec);
        arb_add(m, m, rand + iter, prec);
    }

    arb_div_si(m, m, N, prec);

    for (iter = 0; iter < N; iter++)
    {
        arb_sub(tmp, rand + iter, m, prec);
        arb_sqr(tmp, tmp, prec);
        arb_add(s, s, tmp, prec);
    }

    arb_div_si(s, s, N, prec);

    /* one percent deviation */
    arb_set_str(mp, "0.5 +/- 0.005", prec);
    arb_set_str(sp, "0.083333 +/- 0.00083", prec);

    if (!arb_contains(mp, m))
    {
        flint_printf("FAIL: mean\n\n");
        flint_printf("m = "); arb_printd(m, 15); flint_printf("\n\n");
        flint_abort();
    }

    if (!arb_contains(sp, s))
    {
        flint_printf("FAIL: variance\n\n");
        flint_printf("s = "); arb_printd(s, 15); flint_printf("\n\n");
        flint_abort();
    }

    _arb_vec_clear(rand, N);
    arb_clear(m);
    arb_clear(s);
    arb_clear(mp);
    arb_clear(sp);
    arb_clear(tmp);

    TEST_FUNCTION_END(state);
}
