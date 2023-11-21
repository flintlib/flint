/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_div_2expm1_ui, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        ulong n;
        slong prec, acc1, acc2;
        fmpz_t t;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        fmpz_init(t);

        prec = 2 + n_randint(state, 10000);
        arb_randtest(a, state, 1 + n_randint(state, 10000), 10);

        if (n_randint(state, 2))
            n = 1 + (n_randtest(state) % (10 * prec));
        else
            n = n_randtest(state);

        arb_div_2expm1_ui(b, a, n, prec);

        arb_one(c);
        if (n >= (UWORD(1) << (FLINT_BITS-1)))
        {
            arb_mul_2exp_si(c, c, (UWORD(1) << (FLINT_BITS-2)));
            arb_mul_2exp_si(c, c, (UWORD(1) << (FLINT_BITS-2)));
            arb_mul_2exp_si(c, c, n - (UWORD(1) << (FLINT_BITS-1)));
        }
        else
        {
            arb_mul_2exp_si(c, c, n);
        }

        arb_sub_ui(c, c, 1, prec);
        arb_div(c, a, c, prec);

        acc1 = arb_rel_accuracy_bits(a);
        acc2 = arb_rel_accuracy_bits(b);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("n = %wu\n", n);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        if (n > 0 && acc1 > 1 && (acc2 < FLINT_MIN(prec, acc1) - 10) &&
                !(acc1 == -ARF_PREC_EXACT && acc2 == -ARF_PREC_EXACT))
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec=%wd, acc1=%wd, acc2=%wd\n\n", prec, acc1, acc2);
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
