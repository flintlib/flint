/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_log_rising_ui_naive(acb_t res, const acb_t z, ulong r, slong prec)
{
    acb_t t, u;
    slong k;

    acb_init(t);
    acb_init(u);

    for (k = 0; k < r; k++)
    {
        acb_add_ui(t, z, k, prec);
        acb_log(t, t, prec);
        acb_add(u, u, t, prec);
    }

    acb_swap(res, u);

    acb_clear(t);
    acb_clear(u);
}

TEST_FUNCTION_START(acb_hypgeom_log_rising_ui, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z, s1, s2, a, b;
        slong prec, ebits, prec2;
        ulong n;

        prec = 2 + n_randint(state, 100);

        if (n_randint(state, 10) == 0)
            ebits = 100;
        else
            ebits = 10;
        ebits = 2;

        prec2 = prec + 1 + n_randint(state, 30);

        acb_init(z);
        acb_init(s1);
        acb_init(s2);
        acb_init(a);
        acb_init(b);

        acb_randtest(z, state, prec, ebits);
        acb_randtest(s1, state, prec, 10);
        acb_randtest(s2, state, prec, 10);
        n = n_randint(state, 20);

        if (n_randint(state, 20) == 0)
            n = n_randint(state, 200);

        if (n_randint(state, 2))
        {
            acb_hypgeom_log_rising_ui(s1, z, n, prec);
        }
        else
        {
            acb_set(s1, z);
            acb_hypgeom_log_rising_ui(s1, s1, n, prec);
        }

        acb_hypgeom_log_rising_ui_naive(s2, z, n, prec);

        if (!acb_overlaps(s1, s2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("n = %wu\n\n", n);
            flint_printf("z = "); acb_printn(z, 1000, 0); flint_printf("\n\n");
            flint_printf("s1 = "); acb_printn(s1, 1000, 0); flint_printf("\n\n");
            flint_printf("s2 = "); acb_printn(s2, 1000, 0); flint_printf("\n\n");
            acb_sub(s1, s1, s2, prec2);
            flint_printf("s1 - s2 = "); acb_printd(s1, 1000); flint_printf("\n\n");
            flint_abort();
        }

        acb_get_mid(a, z);

        if (n_randint(state, 2))
        {
            arf_set_mag(arb_midref(acb_realref(b)), arb_radref(acb_realref(z)));
            arf_set_mag(arb_midref(acb_imagref(b)), arb_radref(acb_imagref(z)));

            if (n_randint(state, 2))
                acb_neg(b, b);
            if (n_randint(state, 2))
                acb_conj(b, b);

            acb_add(a, a, b, prec);
        }

        acb_hypgeom_log_rising_ui(s2, a, n, prec);

        if (!acb_overlaps(s1, s2))
        {
            flint_printf("FAIL (2)\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("n = %wu\n\n", n);
            flint_printf("z = "); acb_printn(z, 1000, 0); flint_printf("\n\n");
            flint_printf("a = "); acb_printn(a, 1000, 0); flint_printf("\n\n");
            flint_printf("s1 = "); acb_printn(s1, 1000, 0); flint_printf("\n\n");
            flint_printf("s2 = "); acb_printn(s2, 1000, 0); flint_printf("\n\n");
            acb_sub(s1, s1, s2, prec2);
            flint_printf("s1 - s2 = "); acb_printd(s1, 1000); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z);
        acb_clear(s1);
        acb_clear(s2);
        acb_clear(a);
        acb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
