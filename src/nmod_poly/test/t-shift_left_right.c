/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_shift_left_right, state)
{
    int i, result;

    /* Check a << shift >> shift == a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        slong shift = n_randint(state, 100);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        nmod_poly_shift_left(b, a, shift);
        nmod_poly_shift_right(b, b, shift);

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("shift = %wd, a->length = %wd, n = %wu\n",
                shift, a->length, a->mod.n);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* Check a << shift >> shift == a aliasing the other way */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);
        slong shift = n_randint(state, 100);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(c, state, n_randint(state, 100));

        nmod_poly_set(a, c);
        nmod_poly_shift_left(c, c, shift);
        nmod_poly_shift_right(b, c, shift);

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("shift = %wd, c->length = %wd, n = %wu\n",
                shift, c->length, a->mod.n);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
