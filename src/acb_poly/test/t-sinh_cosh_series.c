/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_sinh_cosh_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t a, b, c, d, e;
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

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);
        acb_poly_init(e);

        acb_poly_randtest(a, state, 20, prec, 3);
        acb_poly_randtest(b, state, 20, prec, 10);
        acb_poly_randtest(c, state, 20, prec, 10);
        acb_poly_randtest(d, state, 20, prec, 10);
        acb_poly_randtest(e, state, 20, prec, 10);

        alg1 = n_randint(state, 4);
        alg2 = n_randint(state, 7);

        switch (alg1)
        {
            case 0:
                acb_poly_sinh_cosh_series(b, c, a, n1, prec);
                break;
            case 1:
                acb_poly_sinh_cosh_series_basecase(b, c, a, n1, prec);
                break;
            case 2:
                acb_poly_sinh_cosh_series_exponential(b, c, a, n1, prec);
                break;
            default:
                acb_poly_sinh_series(b, a, n1, prec);
                acb_poly_cosh_series(c, a, n1, prec);
                break;
        }

        switch (alg2)
        {
            case 0:
                acb_poly_sinh_cosh_series(d, e, a, n2, prec);
                break;
            case 1:
                acb_poly_sinh_cosh_series_basecase(d, e, a, n2, prec);
                break;
            case 2:
                acb_poly_sinh_cosh_series_exponential(d, e, a, n2, prec);
                break;
            case 3:
                acb_poly_sinh_series(d, a, n2, prec);
                acb_poly_cosh_series(e, a, n2, prec);
                break;
            case 4:
                acb_poly_set(d, a);
                acb_poly_sinh_cosh_series(d, e, d, n2, prec);
                break;
            case 5:
                acb_poly_set(e, a);
                acb_poly_sinh_cosh_series(d, e, e, n2, prec);
                break;
            default:
                acb_poly_set(d, a);
                acb_poly_sinh_series(d, d, n2, prec);
                acb_poly_set(e, a);
                acb_poly_cosh_series(e, e, n2, prec);
                break;
        }

        acb_poly_truncate(b, FLINT_MIN(n1, n2));
        acb_poly_truncate(c, FLINT_MIN(n1, n2));
        acb_poly_truncate(d, FLINT_MIN(n1, n2));
        acb_poly_truncate(e, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(b, d) || !acb_poly_overlaps(c, e))
        {
            flint_printf("FAIL\n\n");
            flint_printf("alg1 = %d, alg2 = %d, n1 = %wd, n2 = %wd\n\n", alg1, alg2, n1, n2);

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(e, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
        acb_poly_clear(e);
    }

    TEST_FUNCTION_END(state);
}
