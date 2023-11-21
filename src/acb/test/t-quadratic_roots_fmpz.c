/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_quadratic_roots_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t r1, r2, x, y;
        fmpz_t a, b, c;
        slong prec;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        acb_init(r1);
        acb_init(r2);
        acb_init(x);
        acb_init(y);

        prec = 2 + n_randint(state, 1000);

        fmpz_randtest_not_zero(a, state, 1 + n_randint(state, 1000));
        fmpz_randtest(b, state, 1 + n_randint(state, 1000));
        fmpz_randtest(c, state, 1 + n_randint(state, 1000));

        acb_randtest(r1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(r2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        acb_quadratic_roots_fmpz(r1, r2, a, b, c, prec);

        acb_mul(x, r1, r1, prec);
        acb_mul_fmpz(x, x, a, prec);
        acb_addmul_fmpz(x, r1, b, prec);
        acb_add_fmpz(x, x, c, prec);

        acb_mul(y, r2, r2, prec);
        acb_mul_fmpz(y, y, a, prec);
        acb_addmul_fmpz(y, r2, b, prec);
        acb_add_fmpz(y, y, c, prec);

        if (!acb_contains_zero(x) || !acb_contains_zero(y) ||
            acb_rel_accuracy_bits(r1) < prec - 4 ||
            acb_rel_accuracy_bits(r2) < prec - 4)
        {
            flint_printf("FAIL: containment / accuracy\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf("c = "); fmpz_print(c); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);

        acb_clear(r1);
        acb_clear(r2);
        acb_clear(x);
        acb_clear(y);
    }

    TEST_FUNCTION_END(state);
}
