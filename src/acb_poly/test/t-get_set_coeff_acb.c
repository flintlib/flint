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
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_get_set_coeff_acb, state)
{
    int i, j, result;

    for (i = 0; i < 100; i++)
    {
        acb_poly_t a;
        acb_t x1, x2;
        slong coeff, len;

        acb_poly_init(a);
        acb_init(x1);
        acb_init(x2);
        len = n_randint(state, 100) + 1;

        for (j = 0; j < 100; j++)
        {
            acb_randtest(x1, state, 2 + n_randint(state, 200), 10);
            coeff = n_randint(state, len);
            acb_poly_set_coeff_acb(a, coeff, x1);
            acb_poly_get_coeff_acb(x2, a, coeff);

            result = (acb_equal(x1, x2));
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("x1 = "), acb_print(x1), flint_printf("\n");
                flint_printf("x2 = "), acb_print(x2), flint_printf("\n");
                flint_printf("coeff = %wd, length = %wd\n", coeff, len);
                flint_abort();
            }
        }

        acb_clear(x1);
        acb_clear(x2);
        acb_poly_clear(a);
    }

    TEST_FUNCTION_END(state);
}
