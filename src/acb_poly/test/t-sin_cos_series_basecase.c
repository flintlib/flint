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

    flint_printf("sin_cos_series_basecase....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        slong m, n, rbits1, rbits2;
        fmpq_poly_t B;
        acb_poly_t a, b, c, d, e;
        int times_pi;

        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);
        times_pi = n_randint(state, 2);

        fmpq_poly_init(B);
        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);
        acb_poly_init(e);

        acb_poly_randtest(a, state, m, rbits1, 10);

        acb_poly_sin_cos_series_basecase(b, c, a, n, rbits2, times_pi);

        /* Check sin(x)^2 + cos(x)^2 = 1 */
        acb_poly_mullow(d, b, b, n, rbits2);
        acb_poly_mullow(e, c, c, n, rbits2);
        acb_poly_add(d, d, e, rbits2);

        fmpq_poly_one(B);
        if (!acb_poly_contains_fmpq_poly(d, B))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_set(d, a);
        acb_poly_sin_cos_series_basecase(d, c, d, n, rbits2, times_pi);
        if (!acb_poly_equal(b, d))
        {
            flint_printf("FAIL (aliasing 1)\n\n");
            flint_abort();
        }

        acb_poly_set(d, a);
        acb_poly_sin_cos_series_basecase(b, d, d, n, rbits2, times_pi);
        if (!acb_poly_equal(c, d))
        {
            flint_printf("FAIL (aliasing 2)\n\n");
            flint_abort();
        }

        fmpq_poly_clear(B);
        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
        acb_poly_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

