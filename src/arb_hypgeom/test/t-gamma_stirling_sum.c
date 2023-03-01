/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("gamma_stirling_sum....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t z, s1, s2;
        slong prec, N, K;

        prec = 2 + n_randint(state, 800);
        N = n_randint(state, 200);
        K = n_randint(state, 20);

        arb_init(z);
        arb_init(s1);
        arb_init(s2);

        arb_randtest(z, state, prec, 10 + n_randint(state, 200));

        arb_hypgeom_gamma_stirling_sum_horner(s1, z, N, prec);
        arb_hypgeom_gamma_stirling_sum_improved(s2, z, N, K, prec);

        if (!arb_overlaps(s1, s2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("N = %wd, K = %wd, prec = %wd\n\n", N, K, prec);
            flint_printf("z = "); arb_printn(z, 1000, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 1000, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(z);
        arb_clear(s1);
        arb_clear(s2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

