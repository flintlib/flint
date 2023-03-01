/*
    Copyright (C) 2012 Fredrik Johansson

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

    flint_printf("sinh....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
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
        mpfr_sinh(t, t, MPFR_RNDN);

        arb_sinh(b, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_sinh(a, a, prec);

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

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

