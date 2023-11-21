/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_gamma_stirling_sum, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z, s1, s2;
        slong prec, N, K;

        prec = 2 + n_randint(state, 800);
        N = n_randint(state, 200);
        K = n_randint(state, 20);

        acb_init(z);
        acb_init(s1);
        acb_init(s2);

        acb_randtest(z, state, prec, 10 + n_randint(state, 200));

        acb_hypgeom_gamma_stirling_sum_horner(s1, z, N, prec);
        acb_hypgeom_gamma_stirling_sum_improved(s2, z, N, K, prec);

        if (!acb_overlaps(s1, s2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("N = %wd, K = %wd, prec = %wd\n\n", N, K, prec);
            flint_printf("z = "); acb_printn(z, 1000, 0); flint_printf("\n\n");
            flint_printf("s1 = "); acb_printn(s1, 1000, 0); flint_printf("\n\n");
            flint_printf("s2 = "); acb_printn(s2, 1000, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z);
        acb_clear(s1);
        acb_clear(s2);
    }

    TEST_FUNCTION_END(state);
}
