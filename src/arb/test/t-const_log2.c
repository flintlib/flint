/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arb.h"

TEST_FUNCTION_START(arb_const_log2, state)
{
    slong iter;

    for (iter = 0; iter < 250 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t r;
        mpfr_t s;
        slong accuracy, prec;

        prec = 2 + n_randint(state, 1 << n_randint(state, 16));

        arb_init(r);
        mpfr_init2(s, prec + 1000);

        arb_const_log2(r, prec);
        mpfr_const_log2(s, MPFR_RNDN);

        if (!arb_contains_mpfr(r, s))
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
        mpfr_clear(s);
    }

    TEST_FUNCTION_END(state);
}
