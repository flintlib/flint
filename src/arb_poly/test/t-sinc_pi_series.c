/*
    Copyright (C) 2012, 2013 Fredrik Johansson
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_sinc_pi_series, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n1, n2, rbits1, rbits2, rbits3, rbits4;
        arb_poly_t a, b, c, d;
        arb_t pi;

        rbits1 = 2 + n_randint(state, 300);
        rbits2 = 2 + n_randint(state, 300);
        rbits3 = 2 + n_randint(state, 300);
        rbits4 = 2 + n_randint(state, 300);

        m = n_randint(state, 15);
        n1 = n_randint(state, 15);
        n2 = n_randint(state, 15);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);
        arb_init(pi);

        arb_poly_randtest(a, state, m, rbits1, 10);
        arb_poly_randtest(b, state, 10, rbits1, 10);
        arb_poly_randtest(c, state, 10, rbits1, 10);

        arb_poly_sinc_pi_series(b, a, n1, rbits2);
        arb_poly_sinc_pi_series(c, a, n2, rbits3);

        arb_poly_set(d, b);
        arb_poly_truncate(d, FLINT_MIN(n1, n2));
        arb_poly_truncate(c, FLINT_MIN(n1, n2));

        arb_const_pi(pi, rbits4);

        if (!arb_poly_overlaps(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n1 = %wd, n2 = %wd, bits2 = %wd, bits3 = %wd, bits4 = %wd\n", n1, n2, rbits2, rbits3, rbits4);
            flint_printf("a = "); arb_poly_printd(a, 50); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 50); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 50); flint_printf("\n\n");
            flint_abort();
        }

        /* check pi x sinc_pi(x) = sin_pi(x) */
        arb_poly_mullow(c, b, a, n1, rbits2);
        arb_poly_scalar_mul(c, c, pi, rbits2);
        arb_poly_sin_pi_series(d, a, n1, rbits2);

        if (!arb_poly_overlaps(c, d))
        {
            flint_printf("FAIL (functional equation)\n\n");
            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_poly_printd(d, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_sinc_pi_series(a, a, n1, rbits2);

        if (!arb_poly_overlaps(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
        arb_clear(pi);
    }

    TEST_FUNCTION_END(state);
}
