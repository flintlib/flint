/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_const_reciprocal_fibonacci, state)
{
    slong iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        arb_t r, s, t;
        slong accuracy, prec;

        prec = 2 + n_randint(state, 1 << n_randint(state, 14));

        arb_init(r);
        arb_init(s);
        arb_init(t);

        arb_const_reciprocal_fibonacci(r, prec);
        arb_const_reciprocal_fibonacci(s, prec + 100);
        arb_set_str(t, "3.3598856662431775531720113029189271797 +/- 2.82e-38", 128);

        if (!arb_overlaps(r, s) || !arb_overlaps(r, t) || !arb_overlaps(s, t))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        accuracy = arb_rel_accuracy_bits(r);

        if (accuracy < prec - 4)
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("r = "); arb_printd(r, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(r);
        arb_clear(s);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
