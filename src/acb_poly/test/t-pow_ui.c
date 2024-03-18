/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_pow_ui, state)
{
    slong iter;

    /* compare with fmpz_poly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong zbits1, rbits1, rbits2;
        ulong e;
        fmpz_poly_t A, B;
        acb_poly_t a, b;

        zbits1 = 2 + n_randint(state, 100);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        e = n_randint(state, 30);

        fmpz_poly_init(A);
        fmpz_poly_init(B);

        acb_poly_init(a);
        acb_poly_init(b);

        fmpz_poly_randtest(A, state, 1 + n_randint(state, 8), zbits1);
        fmpz_poly_pow(B, A, e);

        acb_poly_set_fmpz_poly(a, A, rbits1);
        acb_poly_pow_ui(b, a, e, rbits2);

        if (!acb_poly_contains_fmpz_poly(b, B))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);
            flint_printf("e = %wu\n", e);

            flint_printf("A = "); fmpz_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpz_poly_print(B); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_pow_ui(a, a, e, rbits2);
        if (!acb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpz_poly_clear(A);
        fmpz_poly_clear(B);

        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
