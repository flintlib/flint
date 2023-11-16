/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_lgamma_series, state)
{
    slong iter;

    /* special accuracy test case */
    {
        acb_poly_t a;
        acb_t c;

        acb_init(c);
        acb_poly_init(a);

        arb_set_str(acb_realref(c), "-20.25", 53);
        arb_set_str(acb_imagref(c), "1e1000", 53);

        acb_poly_set_coeff_acb(a, 0, c);
        acb_poly_set_coeff_si(a, 1, 1);

        acb_poly_lgamma_series(a, a, 3, 53);

        if (acb_rel_accuracy_bits(a->coeffs) < 40 ||
            acb_rel_accuracy_bits(a->coeffs + 1) < 40 ||
            acb_rel_accuracy_bits(a->coeffs + 2) < 40)
        {
            flint_printf("FAIL: accuracy (reflection formula)\n\n");
            acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_clear(c);
    }

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n1, n2, rbits1, rbits2, rbits3;
        acb_poly_t a, b, c, d;

        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n1 = 1 + n_randint(state, 30);
        n2 = 1 + n_randint(state, 30);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        acb_poly_randtest(a, state, m, rbits1, 10);
        acb_poly_randtest(b, state, m, rbits1, 10);
        acb_poly_randtest(c, state, m, rbits1, 10);

        acb_poly_lgamma_series(b, a, n1, rbits2);
        acb_poly_lgamma_series(c, a, n2, rbits3);

        acb_poly_set(d, b);
        acb_poly_truncate(d, FLINT_MIN(n1, n2));
        acb_poly_truncate(c, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n1 = %wd, n2 = %wd, bits2 = %wd, bits3 = %wd\n", n1, n2, rbits2, rbits3);

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* check loggamma(a) + log(a) = loggamma(a+1) */
        acb_poly_log_series(c, a, n1, rbits2);
        acb_poly_add(c, b, c, rbits2);

        acb_poly_set(d, a);
        acb_add_ui(d->coeffs, d->coeffs, 1, rbits2);
        acb_poly_lgamma_series(d, d, n1, rbits2);

        if (!acb_poly_overlaps(c, d))
        {
            flint_printf("FAIL (functional equation)\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_lgamma_series(a, a, n1, rbits2);
        if (!acb_poly_overlaps(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
