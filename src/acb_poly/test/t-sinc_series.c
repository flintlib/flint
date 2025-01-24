/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_sinc_series, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n1, n2, rbits1, rbits2, rbits3;
        acb_poly_t a, b, c, d;

        rbits1 = 2 + n_randint(state, 300);
        rbits2 = 2 + n_randint(state, 300);
        rbits3 = 2 + n_randint(state, 300);

        m = n_randint(state, 15);
        n1 = n_randint(state, 15);
        n2 = n_randint(state, 15);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        acb_poly_randtest(a, state, m, rbits1, 10);
        acb_poly_randtest(b, state, 10, rbits1, 10);
        acb_poly_randtest(c, state, 10, rbits1, 10);

        acb_poly_sinc_series(b, a, n1, rbits2);
        acb_poly_sinc_series(c, a, n2, rbits3);

        acb_poly_set(d, b);
        acb_poly_truncate(d, FLINT_MIN(n1, n2));
        acb_poly_truncate(c, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n1 = %wd, n2 = %wd, bits2 = %wd, bits3 = %wd\n", n1, n2, rbits2, rbits3);
            flint_printf("a = "); acb_poly_printd(a, 50); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 50); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 50); flint_printf("\n\n");
            flint_abort();
        }

        /* check x sinc(x) = sin(x) */
        acb_poly_mullow(c, b, a, n1, rbits2);
        acb_poly_sin_series(d, a, n1, rbits2);

        if (!acb_poly_overlaps(c, d))
        {
            flint_printf("FAIL (functional equation)\n\n");
            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_sinc_series(a, a, n1, rbits2);

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

    /* check intervals containing 0 */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        acb_t z1, z2;
        acb_poly_t f1, f2, g1, g2;
        slong len, prec;

        acb_init(z1);
        acb_init(z2);

        acb_poly_init(f1);
        acb_poly_init(f2);
        acb_poly_init(g1);
        acb_poly_init(g2);

        prec = 2 + n_randint(state, 300);
        len = 1 + n_randint(state, 5);

        acb_randtest(z1, state, 300, 5);
        acb_mul_2exp_si(z1, z1, -(slong) n_randint(state, 100));
        if (n_randint(state, 2))
            arb_zero(acb_realref(z1));

        arb_get_mag(arb_radref(acb_realref(z2)), acb_realref(z1));
        arb_get_mag(arb_radref(acb_imagref(z2)), acb_imagref(z1));

        acb_poly_set_coeff_si(f1, 1, n_randint(state, 2) ? 1 : -1);
        acb_poly_set(f2, f1);
        acb_poly_set_coeff_acb(f1, 0, z1);
        acb_poly_set_coeff_acb(f2, 0, z2);

        acb_poly_sinc_series(g1, f1, len, prec);
        acb_poly_sinc_series(g2, f2, len, prec);

        if (!acb_poly_overlaps(g1, g2))
        {
            flint_printf("FAIL: overlap (0)\n\n");
            flint_printf("iter = %wd\n", iter);
            flint_printf("f1 = "); acb_poly_printd(f1, prec / 3.33); flint_printf("\n\n");
            flint_printf("f2 = "); acb_poly_printd(f2, prec / 3.33); flint_printf("\n\n");
            flint_printf("g1 = "); acb_poly_printd(g1, prec / 3.33); flint_printf("\n\n");
            flint_printf("g2 = "); acb_poly_printd(g2, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z1);
        acb_clear(z2);

        acb_poly_clear(f1);
        acb_poly_clear(f2);
        acb_poly_clear(g1);
        acb_poly_clear(g2);
    }

    TEST_FUNCTION_END(state);
}
