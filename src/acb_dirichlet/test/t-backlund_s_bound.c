/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_backlund_s_bound, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        mag_t u, v;
        slong aprec, bprec;
        slong abits, bbits;

        aprec = 2 + n_randint(state, 1000);
        bprec = 2 + n_randint(state, 1000);
        abits = 2 + n_randint(state, 100);
        bbits = 2 + n_randint(state, 100);

        arb_init(a);
        arb_init(b);
        mag_init(u);
        mag_init(v);

        arb_randtest(a, state, aprec, abits);
        arb_randtest(b, state, bprec, bbits);

        if (arb_is_nonnegative(a) && arb_is_nonnegative(b))
        {
            acb_dirichlet_backlund_s_bound(u, a);
            acb_dirichlet_backlund_s_bound(v, b);

            if ((arb_lt(a, b) && mag_cmp(u, v) > 0) ||
                (arb_gt(a, b) && mag_cmp(u, v) < 0))
            {
                flint_printf("FAIL: increasing on t >= 0\n\n");
                flint_printf("a = "); arb_print(a); flint_printf("\n\n");
                flint_printf("b = "); arb_print(b); flint_printf("\n\n");
                flint_printf("u = "); mag_print(u); flint_printf("\n\n");
                flint_printf("v = "); mag_print(v); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(a);
        arb_clear(b);
        mag_clear(u);
        mag_clear(v);
    }

    TEST_FUNCTION_END(state);
}
