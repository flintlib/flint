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

    flint_printf("atan2....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
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

        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(r, state, b, 1 + n_randint(state, 200));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        fmpq_get_mpfr(u, r, MPFR_RNDN);
        mpfr_atan2(t, u, t, MPFR_RNDN);

        arb_atan2(c, b, a, prec);

        if (!arb_contains_mpfr(c, t))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        fmpq_clear(q);
        fmpq_clear(r);
        mpfr_clear(t);
        mpfr_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
