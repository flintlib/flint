/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_get_set_coeff_arb, state)
{
    int i, j, result;

    for (i = 0; i < 100; i++)
    {
        arb_poly_t a;
        arb_t x1, x2;
        slong coeff, len;

        arb_poly_init(a);
        arb_init(x1);
        arb_init(x2);
        len = n_randint(state, 100) + 1;

        for (j = 0; j < 100; j++)
        {
            arb_randtest(x1, state, 2 + n_randint(state, 200), 10);
            coeff = n_randint(state, len);
            arb_poly_set_coeff_arb(a, coeff, x1);
            arb_poly_get_coeff_arb(x2, a, coeff);

            result = (arb_equal(x1, x2));
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("x1 = "), arb_print(x1), flint_printf("\n");
                flint_printf("x2 = "), arb_print(x2), flint_printf("\n");
                flint_printf("coeff = %wd, length = %wd\n", coeff, len);
                flint_abort();
            }
        }

        arb_clear(x1);
        arb_clear(x2);
        arb_poly_clear(a);
    }

    TEST_FUNCTION_END(state);
}
