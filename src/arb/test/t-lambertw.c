/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_lambertw, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x1, x2, t, w1, w2;
        slong prec1, prec2, ebits;
        int branch;

        arb_init(x1);
        arb_init(x2);
        arb_init(t);
        arb_init(w1);
        arb_init(w2);

        branch = n_randint(state, 2);

        if (n_randint(state, 4) == 0)
        {
            prec1 = 2 + n_randint(state, 3000);
            prec2 = 2 + n_randint(state, 3000);
            ebits = 1 + n_randint(state, 1000);
        }
        else
        {
            prec1 = 2 + n_randint(state, 300);
            prec2 = 2 + n_randint(state, 300);
            ebits = 1 + n_randint(state, 50);
        }

        arb_randtest(x1, state, 1 + n_randint(state, 1000), ebits);
        arb_randtest(x2, state, 1 + n_randint(state, 1000), ebits);
        arb_randtest(t, state, 1 + n_randint(state, 1000), ebits);
        arb_randtest(w1, state, 1 + n_randint(state, 1000), ebits);
        arb_randtest(w2, state, 1 + n_randint(state, 1000), ebits);

        if (n_randint(state, 4) == 0)
        {
            arb_const_e(t, 2 * prec1);
            arb_inv(t, t, 2 * prec1);
            arb_sub(x1, x1, t, 2 * prec1);
        }

        if (n_randint(state, 2))
        {
            arb_set(x2, x1);
        }
        else
        {
            arb_add(x2, x1, t, 2 * prec1);
            arb_sub(x2, x2, t, 2 * prec1);
        }

        arb_lambertw(w1, x1, branch, prec1);

        if (n_randint(state, 2))
        {
            arb_set(w2, x2);
            arb_lambertw(w2, w2, branch, prec1);
        }
        else
        {
            arb_lambertw(w2, x2, branch, prec1);
        }

        if (!arb_overlaps(w1, w2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("branch = %d, prec1 = %wd, prec2 = %wd\n\n", branch, prec1, prec2);
            flint_printf("x1 = "); arb_printd(x1, 50); flint_printf("\n\n");
            flint_printf("x2 = "); arb_printd(x2, 50); flint_printf("\n\n");
            flint_printf("w1 = "); arb_printd(w1, 50); flint_printf("\n\n");
            flint_printf("w2 = "); arb_printd(w2, 50); flint_printf("\n\n");
            flint_abort();
        }

        arb_exp(t, w1, prec1);
        arb_mul(t, t, w1, prec1);

        if (!arb_contains(t, x1))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("branch = %d, prec1 = %wd, prec2 = %wd\n\n", branch, prec1, prec2);
            flint_printf("x1 = "); arb_printd(x1, 50); flint_printf("\n\n");
            flint_printf("x2 = "); arb_printd(x2, 50); flint_printf("\n\n");
            flint_printf("w1 = "); arb_printd(w1, 50); flint_printf("\n\n");
            flint_printf("w2 = "); arb_printd(w2, 50); flint_printf("\n\n");
            flint_printf("t  = "); arb_printd(t, 50); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x1);
        arb_clear(x2);
        arb_clear(t);
        arb_clear(w1);
        arb_clear(w2);
    }

    TEST_FUNCTION_END(state);
}
