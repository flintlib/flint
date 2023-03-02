/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("acos_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A;
        arb_poly_t a, b, c, d;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);

        fmpq_poly_init(A);
        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        fmpq_poly_randtest(A, state, m, qbits);
        arb_poly_set_fmpq_poly(a, A, rbits1);

        arb_poly_randtest(b, state, 1 + n_randint(state, 30), rbits1, 5);

        arb_poly_acos_series(b, a, n, rbits2);

        /* Check cos(acos(x)) = x */
        arb_poly_sin_cos_series_basecase(d, c, b, n, rbits2, 0);

        fmpq_poly_truncate(A, n);
        if (!arb_poly_contains_fmpq_poly(c, A))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_acos_series(a, a, n, rbits2);
        if (!arb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpq_poly_clear(A);
        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

