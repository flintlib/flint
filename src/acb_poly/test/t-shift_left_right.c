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

TEST_FUNCTION_START(acb_poly_shift_left_right, state)
{
    int i, result;

    /* Check aliasing of a and b for left shift */
    for (i = 0; i < 1000; i++)
    {
        acb_poly_t a, b;
        slong shift = n_randint(state, 100);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        acb_poly_shift_left(b, a, shift);
        acb_poly_shift_left(a, a, shift);

        result = (acb_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            acb_poly_printd(a, 10), flint_printf("\n\n");
            acb_poly_printd(b, 10), flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    /* Check aliasing of a and b for right shift */
    for (i = 0; i < 1000; i++)
    {
        acb_poly_t a, b;
        slong shift = n_randint(state, 100);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        acb_poly_shift_right(b, a, shift);
        acb_poly_shift_right(a, a, shift);

        result = (acb_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            acb_poly_printd(a, 10), flint_printf("\n\n");
            acb_poly_printd(b, 10), flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    /* Check shift left then right does nothing */
    for (i = 0; i < 1000; i++)
    {
        acb_poly_t a, b, c;
        slong shift = n_randint(state, 100);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        acb_poly_shift_left(b, a, shift);
        acb_poly_shift_right(c, b, shift);

        result = (acb_poly_equal(c, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            acb_poly_printd(a, 10), flint_printf("\n\n");
            acb_poly_printd(b, 10), flint_printf("\n\n");
            acb_poly_printd(c, 10), flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
