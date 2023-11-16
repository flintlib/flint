/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_sinh_cosh_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_poly_t a, b, c, d, e;
        slong n1, n2, prec;
        int alg1, alg2;

        prec = 2 + n_randint(state, 200);

        if (n_randint(state, 40) == 0)
        {
            n1 = n_randint(state, 300);
            n2 = n_randint(state, 300);
        }
        else
        {
            n1 = n_randint(state, 30);
            n2 = n_randint(state, 30);
        }

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);
        arb_poly_init(e);

        arb_poly_randtest(a, state, 20, prec, 10);
        arb_poly_randtest(b, state, 20, prec, 10);
        arb_poly_randtest(c, state, 20, prec, 10);
        arb_poly_randtest(d, state, 20, prec, 10);
        arb_poly_randtest(e, state, 20, prec, 10);

        alg1 = n_randint(state, 4);
        alg2 = n_randint(state, 7);

        switch (alg1)
        {
            case 0:
                arb_poly_sinh_cosh_series(b, c, a, n1, prec);
                break;
            case 1:
                arb_poly_sinh_cosh_series_basecase(b, c, a, n1, prec);
                break;
            case 2:
                arb_poly_sinh_cosh_series_exponential(b, c, a, n1, prec);
                break;
            default:
                arb_poly_sinh_series(b, a, n1, prec);
                arb_poly_cosh_series(c, a, n1, prec);
                break;
        }

        switch (alg2)
        {
            case 0:
                arb_poly_sinh_cosh_series(d, e, a, n2, prec);
                break;
            case 1:
                arb_poly_sinh_cosh_series_basecase(d, e, a, n2, prec);
                break;
            case 2:
                arb_poly_sinh_cosh_series_exponential(d, e, a, n2, prec);
                break;
            case 3:
                arb_poly_sinh_series(d, a, n2, prec);
                arb_poly_cosh_series(e, a, n2, prec);
                break;
            case 4:
                arb_poly_set(d, a);
                arb_poly_sinh_cosh_series(d, e, d, n2, prec);
                break;
            case 5:
                arb_poly_set(e, a);
                arb_poly_sinh_cosh_series(d, e, e, n2, prec);
                break;
            default:
                arb_poly_set(d, a);
                arb_poly_sinh_series(d, d, n2, prec);
                arb_poly_set(e, a);
                arb_poly_cosh_series(e, e, n2, prec);
                break;
        }

        arb_poly_truncate(b, FLINT_MIN(n1, n2));
        arb_poly_truncate(c, FLINT_MIN(n1, n2));
        arb_poly_truncate(d, FLINT_MIN(n1, n2));
        arb_poly_truncate(e, FLINT_MIN(n1, n2));

        if (!arb_poly_overlaps(b, d) || !arb_poly_overlaps(c, e))
        {
            flint_printf("FAIL\n\n");
            flint_printf("alg1 = %d, alg2 = %d, n1 = %wd, n2 = %wd\n\n", alg1, alg2, n1, n2);

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_poly_printd(d, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_poly_printd(e, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
        arb_poly_clear(e);
    }

    TEST_FUNCTION_END(state);
}
