/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

#define N 10000

/* same as t-urandom in arb/, but ignore variance */
int main(void)
{
    slong iter;
    slong prec;
    flint_rand_t state;
    acb_ptr rand;
    acb_t m; /* mean */
    acb_t mp;
    arb_t tol;

    flint_printf("urandom....");
    fflush(stdout);

    flint_randinit(state);
    acb_init(m);
    acb_init(mp);
    arb_init(tol);
    rand = _acb_vec_init(N);
    prec = 299;

    for (iter = 0; iter < N; iter++)
    {
        acb_urandom(rand + iter, state, prec);
        acb_add(m, m, rand + iter, prec);
    }

    acb_div_si(m, m, N, prec);

    /* one percent deviation */
    arb_set_str(tol, "0 +/- 0.01", prec);
    acb_set_arb_arb(mp, tol, tol);

    if (!acb_contains(mp, m))
    {
        flint_printf("FAIL: mean\n\n");
        flint_printf("m = "); acb_printd(m, 15); flint_printf("\n\n");
        flint_abort();
    }

    _acb_vec_clear(rand, N);
    acb_clear(m);
    acb_clear(mp);
    arb_clear(tol);
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
