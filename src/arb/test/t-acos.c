/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("acos....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t q;
        mpfr_t t;
        slong prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        fmpq_init(q);
        mpfr_init2(t, prec + 100);

        arb_randtest(a, state, 1 + n_randint(state, 200), 3);
        arb_randtest(b, state, 1 + n_randint(state, 200), 3);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_acos(t, t, MPFR_RNDN);

        arb_acos(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_printd(a, 100); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 100); flint_printf("\n\n");
            flint_abort();
        }

        arb_acos(a, a, prec);

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

    /* check large arguments */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_precise(a, state, 1 + n_randint(state, 1000), 100);

        arb_acos(b, a, prec1);
        arb_acos(c, a, prec2);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        /* check sin(asin(x)) = x */
        arb_cos(c, b, prec1);

        if (!arb_contains(c, a))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
