/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_sub, state)
{
    int i, result;

    /* Check a - b = a + neg(b) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        nmod_poly_sub(c, a, b);
        nmod_poly_neg(b, b);
        nmod_poly_add(d, a, b);

        result = (nmod_poly_equal(d, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            nmod_poly_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        nmod_poly_sub(c, a, b);
        nmod_poly_sub(a, a, b);

        result = (nmod_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        nmod_poly_sub(c, a, b);
        nmod_poly_sub(b, a, b);

        result = (nmod_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
