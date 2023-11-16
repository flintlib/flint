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
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_agm, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        fmpq_t q, r;
        mpfr_t t, u;
        slong prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        fmpq_init(q);
        fmpq_init(r);
        mpfr_init2(t, prec + 100);
        mpfr_init2(u, prec + 100);

        arb_randtest(a, state, 1 + n_randint(state, 200), 3);
        arb_randtest(b, state, 1 + n_randint(state, 200), 3);
        arb_randtest(c, state, 1 + n_randint(state, 200), 3);

        arb_agm(c, a, b, prec);

        if (arb_equal(a, b))
        {
            if (!arb_contains(c, a))
            {
                flint_printf("FAIL: containment (identity)\n\n");
                flint_printf("a = "); arb_print(a); flint_printf("\n\n");
                flint_printf("b = "); arb_print(b); flint_printf("\n\n");
                flint_printf("c = "); arb_print(c); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));
            arb_get_rand_fmpq(r, state, b, 1 + n_randint(state, 200));
            fmpq_get_mpfr(t, q, MPFR_RNDN);
            fmpq_get_mpfr(u, r, MPFR_RNDN);
            mpfr_agm(t, t, u, MPFR_RNDN);

            if (!arb_contains_mpfr(c, t))
            {
                flint_printf("FAIL: containment\n\n");
                flint_printf("a = "); arb_print(a); flint_printf("\n\n");
                flint_printf("b = "); arb_print(b); flint_printf("\n\n");
                flint_printf("c = "); arb_print(c); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        fmpq_clear(q);
        fmpq_clear(r);
        mpfr_clear(t);
        mpfr_clear(u);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x1, x2, y1, y2, r1, r2;
        slong prec1, prec2;

        arb_init(x1);
        arb_init(x2);
        arb_init(y1);
        arb_init(y2);
        arb_init(r1);
        arb_init(r2);

        arb_randtest_special(x1, state, 1 + n_randint(state, 200), 100);
        arb_randtest_special(y1, state, 1 + n_randint(state, 200), 100);

        if (n_randint(state, 2))
        {
            arb_randtest_special(r1, state, 1 + n_randint(state, 200), 100);
            arb_randtest_special(r2, state, 1 + n_randint(state, 200), 100);
        }
        else
        {
            arb_zero(r1);
            arb_zero(r2);
        }

        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);

        arb_add(x2, x1, r1, prec2);
        arb_sub(x2, x2, r1, prec2);

        arb_add(y2, y1, r2, prec2);
        arb_sub(y2, y2, r2, prec2);

        arb_agm(r1, x1, y1, prec1);
        arb_agm(r2, x2, y2, prec2);

        if (!arb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x1 = "); arb_printn(x1, 30, 0); flint_printf("\n\n");
            flint_printf("x2 = "); arb_printn(x2, 30, 0); flint_printf("\n\n");
            flint_printf("y1 = "); arb_printn(y1, 30, 0); flint_printf("\n\n");
            flint_printf("y2 = "); arb_printn(y2, 30, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 30, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 30, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x1);
        arb_clear(x2);
        arb_clear(y1);
        arb_clear(y2);
        arb_clear(r1);
        arb_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
