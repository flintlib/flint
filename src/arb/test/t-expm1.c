/*
    Copyright (C) 2012, 2013 Fredrik Johansson

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

TEST_FUNCTION_START(arb_expm1, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        slong prec0, prec;

        prec0 = 400;
        if (iter % 100 == 0)
            prec0 = 10000;

        prec = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);

        arb_randtest(a, state, 1 + n_randint(state, prec0), 4);
        arb_randtest(b, state, 1 + n_randint(state, prec0), 4);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, prec0));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_expm1(t, t, MPFR_RNDN);

        arb_expm1(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("iter = %wd, prec = %wd\n\n", iter, prec);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_expm1(a, a, prec);

        if (!arb_equal(a, b))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("iter = %wd, prec = %wd\n\n", iter, prec);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        fmpq_clear(q);
        mpfr_clear(t);
    }

    /* check large arguments */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        slong prec0, prec1, prec2, acc1, acc2;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 1000;

        prec1 = 2 + n_randint(state, prec0);
        prec2 = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);

        arb_expm1(b, a, prec1);
        arb_set(c, a);  /* also tests aliasing */
        arb_expm1(c, c, prec2);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);

        /* compare with exp */
        arb_expm1(c, a, prec1);

        arb_exp(d, a, prec1);
        arb_sub_ui(d, d, 1, prec1);

        acc1 = arb_rel_accuracy_bits(c);
        acc2 = arb_rel_accuracy_bits(d);

        if (!arb_overlaps(c, d) || ((acc2 > 0) && acc1 < FLINT_MIN(acc2, prec1) - 3))
        {
            flint_printf("FAIL: comparison with exp\n\n");
            flint_printf("prec = %wd, acc1 = %wd, acc2 = %wd\n\n", prec1, acc1, acc2);
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
