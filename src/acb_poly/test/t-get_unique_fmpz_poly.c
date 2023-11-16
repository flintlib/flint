/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_get_unique_fmpz_poly, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec, c;
        fmpz_poly_t A, B, C;
        acb_poly_t a, b;

        fmpz_poly_init(A);
        fmpz_poly_init(B);
        fmpz_poly_init(C);
        acb_poly_init(a);
        acb_poly_init(b);

        fmpz_poly_randtest(A, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000));
        fmpz_poly_randtest(B, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000));
        fmpz_poly_randtest(C, state, 1 + n_randint(state, 10), 1 + n_randint(state, 1000));
        c = 1 + n_randint(state, 1000);

        prec = 2 + n_randint(state, 100);

        for ( ; ; )
        {
            acb_poly_set_fmpz_poly(a, A, prec);
            acb_poly_set2_fmpz_poly(b, B, C, prec);
            acb_poly_scalar_mul_2exp_si(b, b, -c);
            acb_poly_add(a, a, b, prec);
            acb_poly_sub(a, a, b, prec);

            if (acb_poly_get_unique_fmpz_poly(B, a))
            {
                if (!fmpz_poly_equal(A, B))
                {
                    flint_printf("FAIL\n\n");
                    flint_printf("A = "); fmpz_poly_print(A); flint_printf("\n\n");
                    flint_printf("B = "); fmpz_poly_print(B); flint_printf("\n\n");
                    flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
                    flint_abort();
                }

                break;
            }
            else
            {
                prec *= 2;
            }
        }

        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
        fmpz_poly_clear(C);
        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
