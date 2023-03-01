/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("rsqrt_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 5000 * arb_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A;
        acb_poly_t a, b, c;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);

        fmpq_poly_init(A);
        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);

        fmpq_poly_randtest_not_zero(A, state, m, qbits);
        acb_poly_set_fmpq_poly(a, A, rbits1);

        acb_poly_randtest(b, state, 1 + n_randint(state, 30), rbits1, 5);

        acb_poly_rsqrt_series(b, a, n, rbits2);

        /* Check 1/((1/sqrt(a))^2) = a */
        acb_poly_mullow(c, b, b, n, rbits2);
        acb_poly_inv_series(c, c, n, rbits2);

        fmpq_poly_truncate(A, n);
        if (!acb_poly_contains_fmpq_poly(c, A))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_rsqrt_series(a, a, n, rbits2);
        if (!acb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpq_poly_clear(A);
        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

