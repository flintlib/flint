/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_lambertw_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t x1, x2, t, w1, w2;
        slong prec1, prec2, n1, n2, ebits;
        fmpz_t k;

        acb_poly_init(x1);
        acb_poly_init(x2);
        acb_poly_init(t);
        acb_poly_init(w1);
        acb_poly_init(w2);
        fmpz_init(k);

        n1 = n_randint(state, 15);
        n2 = n_randint(state, 15);

        if (n_randint(state, 4) == 0)
        {
            prec1 = 2 + n_randint(state, 3000);
            prec2 = 2 + n_randint(state, 3000);
            ebits = 1 + n_randint(state, 1000);
        }
        else
        {
            prec1 = 2 + n_randint(state, 300);
            prec2 = 2 + n_randint(state, 300);
            ebits = 1 + n_randint(state, 50);
        }

        fmpz_randtest(k, state, ebits);
        acb_poly_randtest(x1, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000), ebits);
        acb_poly_randtest(x2, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000), ebits);
        acb_poly_randtest(t, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000), ebits);
        acb_poly_randtest(w1, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000), ebits);
        acb_poly_randtest(w2, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000), ebits);

        if (n_randint(state, 4) == 0)
        {
            acb_t u;
            acb_init(u);
            arb_const_e(acb_realref(u), 2 * prec1);
            acb_inv(u, u, 2 * prec1);
            acb_poly_set_acb(t, u);
            acb_poly_sub(x1, x1, t, 2 * prec1);
            acb_clear(u);
        }

        if (n_randint(state, 2))
        {
            acb_poly_set(x2, x1);
        }
        else
        {
            acb_poly_add(x2, x1, t, 2 * prec1);
            acb_poly_sub(x2, x2, t, 2 * prec1);
        }

        acb_poly_lambertw_series(w1, x1, k, 0, n1, prec1);

        if (n_randint(state, 2))
        {
            acb_poly_set(w2, x2);
            acb_poly_lambertw_series(w2, w2, k, 0, n2, prec1);
        }
        else
        {
            acb_poly_lambertw_series(w2, x2, k, 0, n2, prec1);
        }

        acb_poly_truncate(w1, FLINT_MIN(n1, n2));
        acb_poly_truncate(w2, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(w1, w2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter %wd, branch = ", iter); fmpz_print(k);
            flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
            flint_printf("x1 = "); acb_poly_printd(x1, 50); flint_printf("\n\n");
            flint_printf("x2 = "); acb_poly_printd(x2, 50); flint_printf("\n\n");
            flint_printf("w1 = "); acb_poly_printd(w1, 50); flint_printf("\n\n");
            flint_printf("w2 = "); acb_poly_printd(w2, 50); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_exp_series(t, w1, FLINT_MIN(n1, n2), prec1);
        acb_poly_mullow(t, t, w1, FLINT_MIN(n1, n2), prec1);
        acb_poly_truncate(x1, FLINT_MIN(n1, n2));

        if (!acb_poly_contains(t, x1))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("iter %wd, branch = ", iter); fmpz_print(k);
            flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
            flint_printf("x1 = "); acb_poly_printd(x1, 50); flint_printf("\n\n");
            flint_printf("x2 = "); acb_poly_printd(x2, 50); flint_printf("\n\n");
            flint_printf("w1 = "); acb_poly_printd(w1, 50); flint_printf("\n\n");
            flint_printf("w2 = "); acb_poly_printd(w2, 50); flint_printf("\n\n");
            flint_printf("t  = "); acb_poly_printd(t, 50); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(x1);
        acb_poly_clear(x2);
        acb_poly_clear(t);
        acb_poly_clear(w1);
        acb_poly_clear(w2);
        fmpz_clear(k);
    }

    TEST_FUNCTION_END(state);
}
