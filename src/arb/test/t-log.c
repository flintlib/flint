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

TEST_FUNCTION_START(arb_log, state)
{
    slong iter;

    /* compare with mpfr */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        slong prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        } while (arb_contains_nonpositive(a));

        arb_randtest(b, state, 1 + n_randint(state, 200), 10);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);

        /* todo: estimate cancellation precisely */
        if (mpfr_cmp_d(t, 1 - 1e-10) > 0 && mpfr_cmp_d(t, 1 + 1e-10) < 0)
        {
            mpfr_set_prec(t, prec + 1000);
            fmpq_get_mpfr(t, q, MPFR_RNDN);
        }

        mpfr_log(t, t, MPFR_RNDN);

        arb_log(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_log(a, a, prec);

        if (!arb_equal(a, b))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        fmpq_clear(q);
        mpfr_clear(t);
    }

    /* compare with mpfr (higher precision) */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        slong prec = 2 + n_randint(state, 6000);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 6000), 10);
        } while (arb_contains_nonpositive(a));

        if (n_randint(state, 10) == 0)
            arb_set_ui(a, n_randint(state, 100));

        arb_randtest(b, state, 1 + n_randint(state, 6000), 10);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);

        /* todo: estimate cancellation precisely */
        if (mpfr_cmp_d(t, 1 - 1e-10) > 0 && mpfr_cmp_d(t, 1 + 1e-10) < 0)
        {
            mpfr_set_prec(t, prec + 10000);
            fmpq_get_mpfr(t, q, MPFR_RNDN);
        }

        mpfr_log(t, t, MPFR_RNDN);

        arb_log(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_log(a, a, prec);

        if (!arb_equal(a, b))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        fmpq_clear(q);
        mpfr_clear(t);
    }

    /* test large numbers */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, ab, lab, la, lb, lalb;
        slong prec = 2 + n_randint(state, 6000);

        arb_init(a);
        arb_init(b);
        arb_init(ab);
        arb_init(lab);
        arb_init(la);
        arb_init(lb);
        arb_init(lalb);

        arb_randtest(a, state, 1 + n_randint(state, 400), 400);
        arb_randtest(b, state, 1 + n_randint(state, 400), 400);

        arb_log(la, a, prec);
        arb_log(lb, b, prec);
        arb_mul(ab, a, b, prec);
        arb_log(lab, ab, prec);
        arb_add(lalb, la, lb, prec);

        if (!arb_overlaps(lab, lalb))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("la = "); arb_print(la); flint_printf("\n\n");
            flint_printf("lb = "); arb_print(lb); flint_printf("\n\n");
            flint_printf("ab = "); arb_print(ab); flint_printf("\n\n");
            flint_printf("lab = "); arb_print(lab); flint_printf("\n\n");
            flint_printf("lalb = "); arb_print(lalb); flint_printf("\n\n");
            flint_abort();
        }

        arb_log(a, a, prec);
        if (!arb_overlaps(a, la))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(ab);
        arb_clear(lab);
        arb_clear(la);
        arb_clear(lb);
        arb_clear(lalb);
    }

    TEST_FUNCTION_END(state);
}
