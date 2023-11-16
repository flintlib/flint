/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_pow_fmpq, state)
{
    slong iter;

    /* check large arguments */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        fmpq_t e1, e2, e3;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        fmpq_init(e1);
        fmpq_init(e2);
        fmpq_init(e3);

        arb_randtest_precise(a, state, 1 + n_randint(state, 1000), 200);
        arb_randtest_precise(b, state, 1 + n_randint(state, 1000), 200);
        fmpq_randtest(e1, state, 200);
        fmpq_randtest(e2, state, 200);

        arb_pow_fmpq(b, a, e1, prec1);
        arb_pow_fmpq(c, a, e1, prec2);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("e1 = "); fmpq_print(e1); flint_printf("\n\n");
            flint_abort();
        }

        /* check a^(e1+e2) = a^e1*a^e2 */
        arb_pow_fmpq(c, a, e2, prec1);
        arb_mul(d, b, c, prec1);
        fmpq_add(e3, e1, e2);
        arb_pow_fmpq(c, a, e3, prec1);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_printf("e1 = "); fmpq_print(e1); flint_printf("\n\n");
            flint_printf("e2 = "); fmpq_print(e2); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        fmpq_clear(e1);
        fmpq_clear(e2);
        fmpq_clear(e3);
    }

    TEST_FUNCTION_END(state);
}
