/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_roots_fmpz_poly, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t poly, f;
        qqbar_ptr r;
        acb_t x, y;
        slong i, d, prec;

        fmpz_poly_init(poly);
        fmpz_poly_init(f);

        fmpz_poly_randtest(poly, state, 10, 100);
        if (n_randint(state, 2))
        {
            fmpz_poly_randtest_not_zero(f, state, 5, 10);
            fmpz_poly_pow(f, f, 1 + n_randint(state, 3));
            fmpz_poly_mul(poly, poly, f);
        }

        d = fmpz_poly_degree(poly);
        d = FLINT_MAX(d, 0);

        acb_init(x);
        acb_init(y);
        r = _qqbar_vec_init(d);

        qqbar_roots_fmpz_poly(r, poly, 0);

        for (i = 0; i < d; i++)
        {
            prec = 2 * QQBAR_DEFAULT_PREC;

            qqbar_get_acb(x, r + i, prec);
            arb_fmpz_poly_evaluate_acb(y, poly, x, prec);

            if (!acb_contains_zero(y))
            {
                flint_printf("FAIL!\n");
                flint_printf("i = %wd\n\n", i);
                flint_printf("poly = "); fmpz_poly_print(poly); flint_printf("\n\n");
                flint_printf("r = "); qqbar_print(r + i); flint_printf("\n\n");
                flint_printf("x = "); acb_printn(x, 100, 0); flint_printf("\n\n");
                flint_printf("y = "); acb_printn(y, 100, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        fmpz_poly_clear(poly);
        fmpz_poly_clear(f);
        _qqbar_vec_clear(r, d);
        acb_clear(x);
        acb_clear(y);
    }

    TEST_FUNCTION_END(state);
}
