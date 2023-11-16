/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_exp_invexp, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        slong prec, acc1, acc2;

        if (iter % 10 == 0)
            prec = 10000;
        else
            prec = 1000;

        prec = 2 + n_randint(state, prec);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, prec), 1 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, prec), 1 + n_randint(state, 100));
        arb_randtest_special(c, state, 1 + n_randint(state, prec), 1 + n_randint(state, 100));

        if (n_randint(state, 2))
        {
            arb_exp_invexp(b, c, a, prec);
        }
        else if (n_randint(state, 2))
        {
            arb_set(b, a);
            arb_exp_invexp(b, c, b, prec);
        }
        else
        {
            arb_set(c, a);
            arb_exp_invexp(b, c, c, prec);
        }

        arb_exp(d, a, prec);

        acc1 = arb_rel_accuracy_bits(b);
        acc2 = arb_rel_accuracy_bits(d);

        if (!arb_overlaps(b, d) || ((acc1 > 0 || acc2 > 0) && acc1 < FLINT_MIN(acc2, prec) - 3))
        {
            flint_printf("FAIL: overlap 1\n\n");
            flint_printf("prec = %wd, acc1 = %wd, acc2 = %wd\n\n", prec, acc1, acc2);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_neg(d, a);
        arb_exp(d, d, prec);

        acc1 = arb_rel_accuracy_bits(c);
        acc2 = arb_rel_accuracy_bits(d);

        if (!arb_overlaps(c, d) || ((acc1 > 0 || acc2 > 0) && acc1 < FLINT_MIN(acc2, prec) - 3))
        {
            flint_printf("FAIL: overlap 2\n\n");
            flint_printf("prec = %wd, acc1 = %wd, acc2 = %wd\n\n", prec, acc1, acc2);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
