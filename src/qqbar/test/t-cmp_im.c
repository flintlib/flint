/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_cmp_im, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, xr, yr, t, u;
        int c1, c2;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(xr);
        qqbar_init(yr);
        qqbar_init(t);
        qqbar_init(u);

        qqbar_randtest(x, state, 3, 100);
        if (n_randint(state, 2))
            qqbar_set(y, x);
        else
            qqbar_randtest(y, state, 3, 10);
        if (n_randint(state, 2))
            qqbar_conj(y, y);

        qqbar_randtest_real(t, state, 1, 100);
        qqbar_i(u);
        qqbar_mul(t, t, u);
        qqbar_add(x, x, t);

        qqbar_im(xr, x);
        qqbar_im(yr, y);

        c1 = qqbar_cmp_im(x, y);
        c2 = qqbar_cmp_re(xr, yr);

        if (c1 != c2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("%d\n\n", c1);
            flint_printf("%d\n\n", c2);
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(xr);
        qqbar_clear(yr);
        qqbar_clear(t);
        qqbar_clear(u);
    }

    /* Some branches are difficult to reach with random test cases. */
    {
        qqbar_t x, y, z, i;
        mag_t eps;
        int which, ans, want;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(i);
        mag_init(eps);

        qqbar_i(i);

        for (which = 0; which <= 8; which++)
        {
            if (which <= 4)
            {
                qqbar_set_si(x, 5);
                qqbar_mul_si(z, i, 3);
                qqbar_add(x, x, z);
                qqbar_set_si(y, 6);
                qqbar_mul_si(z, i, 3);
                qqbar_add(y, y, z);

                if (which == 0)
                {
                    qqbar_set_d(z, 1e-30);
                    qqbar_mul(z, z, i);
                    qqbar_add(x, x, z);
                    want = 1;
                }
                else if (which == 1)
                {
                    qqbar_set_d(z, -1e-30);
                    qqbar_mul(z, z, i);
                    qqbar_add(x, x, z);
                    want = -1;
                }
                else if (which == 2)
                {
                    qqbar_set_d(z, -1e-30);
                    qqbar_mul(z, z, i);
                    qqbar_add(x, x, z);
                    mag_set_d(eps, 1e-15);
                    arb_add_error_mag(acb_imagref(QQBAR_ENCLOSURE(y)), eps);
                    want = -1;
                }
                else if (which == 3)
                {
                    qqbar_set_d(z, 1e-30);
                    qqbar_mul(z, z, i);
                    qqbar_add(x, x, z);
                    mag_set_d(eps, 1e-15);
                    arb_add_error_mag(acb_imagref(QQBAR_ENCLOSURE(y)), eps);
                    want = 1;
                }
                else if (which == 4)
                {
                    mag_set_d(eps, 1e-10);
                    arb_add_error_mag(acb_imagref(QQBAR_ENCLOSURE(x)), eps);
                    mag_set_d(eps, 1e-15);
                    arb_add_error_mag(acb_imagref(QQBAR_ENCLOSURE(y)), eps);
                    want = 0;
                }
            }
            else if (which == 5)
            {
                qqbar_set_si(x, 5);
                qqbar_set_si(y, 6);
                mag_set_d(eps, 1e-10);
                arb_add_error_mag(acb_imagref(QQBAR_ENCLOSURE(y)), eps);
                want = 0;
            }
            else if (which == 6)
            {
                qqbar_sqrt_ui(x, 2);
                qqbar_set_si(y, 1);
                qqbar_add(y, y, i);
                mag_set_d(eps, 1.1);
                arb_add_error_mag(acb_imagref(QQBAR_ENCLOSURE(y)), eps);
                want = -1;
            }
            else if (which == 7)
            {
                qqbar_sqrt_ui(x, 2);
                qqbar_set_si(y, 1);
                qqbar_sub(y, y, i);
                mag_set_d(eps, 1.1);
                arb_add_error_mag(acb_imagref(QQBAR_ENCLOSURE(y)), eps);
                want = 1;
            }
            else if (which == 8)
            {
                qqbar_sqrt_ui(x, 2);
                qqbar_sqrt_ui(y, 3);
                mag_set_d(eps, 0.5);
                arb_add_error_mag(acb_imagref(QQBAR_ENCLOSURE(y)), eps);
                want = 0;
            }

            ans = qqbar_cmp_im(x, y);
     
            if (ans != want)
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
                flint_printf("%d\n\n", ans);
                flint_printf("%d\n\n", want);
                flint_abort();
            }
       }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(i);
        mag_clear(eps);
    }

    TEST_FUNCTION_END(state);
}
