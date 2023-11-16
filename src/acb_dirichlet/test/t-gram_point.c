/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_gram_point, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x1, x2, v1, v2, pin;
        acb_t t;
        fmpz_t n;
        slong prec1, prec2;

        arb_init(x1);
        arb_init(x2);
        arb_init(v1);
        arb_init(v2);
        arb_init(pin);
        acb_init(t);
        fmpz_init(n);

        fmpz_randtest_unsigned(n, state, 500);
        fmpz_sub_ui(n, n, 1);
        prec1 = 2 + n_randtest(state) % 500;
        prec2 = 2 + n_randtest(state) % 2000;

        acb_dirichlet_gram_point(x1, n, NULL, NULL, prec1);
        acb_dirichlet_gram_point(x2, n, NULL, NULL, prec2);

        arb_const_pi(pin, FLINT_MAX(prec1, prec2) + 20);
        arb_mul_fmpz(pin, pin, n, FLINT_MAX(prec1, prec2) + 20);

        acb_set_arb(t, x1);
        acb_dirichlet_hardy_theta(t, t, NULL, NULL, 1, prec1 + 20);
        arb_set(v1, acb_realref(t));

        acb_set_arb(t, x2);
        acb_dirichlet_hardy_theta(t, t, NULL, NULL, 1, prec2 + 20);
        arb_set(v2, acb_realref(t));

        if (!arb_overlaps(x1, x2) || !arb_contains(v1, pin) || !arb_contains(v2, pin))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = "); fmpz_print(n);
            flint_printf("   prec1 = %wd  prec2 = %wd\n\n", prec1, prec2);
            flint_printf("x1 = "); arb_printn(x1, 100, 0); flint_printf("\n\n");
            flint_printf("x2 = "); arb_printn(x2, 100, 0); flint_printf("\n\n");
            flint_printf("v1 = "); arb_printn(v1, 100, 0); flint_printf("\n\n");
            flint_printf("v2 = "); arb_printn(v2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(x1) < prec1 - 3 || arb_rel_accuracy_bits(x2) < prec2 - 3)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("n = "); fmpz_print(n);
            flint_printf("   prec1 = %wd  prec2 = %wd\n\n", prec1, prec2);
            flint_printf("acc(x1) = %wd, acc(x2) = %wd\n\n", arb_rel_accuracy_bits(x1), arb_rel_accuracy_bits(x2));
            flint_printf("x1 = "); arb_printn(x1, 100, 0); flint_printf("\n\n");
            flint_printf("x2 = "); arb_printn(x2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x1);
        arb_clear(x2);
        arb_clear(v1);
        arb_clear(v2);
        arb_clear(pin);
        acb_clear(t);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
