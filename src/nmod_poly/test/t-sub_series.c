/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2009, 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_sub_series, state)
{
    int i, result;

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p;
        nmod_poly_t a, b, c;
        slong n;

        p = n_randtest(state);
        if (p == 0)
	   p++;

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(c, p);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        n = n_randint(state, 100);

        nmod_poly_sub_series(c, a, b, n);
        nmod_poly_sub_series(a, a, b, n);

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
        ulong p;
        nmod_poly_t a, b, c;
        slong n;

        p = n_randtest(state);
        if (p == 0)
	   p++;

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(c, p);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        n = n_randint(state, 100);

        nmod_poly_sub_series(c, a, b, n);
        nmod_poly_sub_series(b, a, b, n);

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

    /* Check truncate(a + b, n) = add_series(a, b, n) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p;
        nmod_poly_t a, b, c, d;
        slong n;

        p = n_randtest(state);
        if (p == 0)
           p++;

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(c, p);
        nmod_poly_init(d, p);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        n = n_randint(state, 100);

        nmod_poly_sub(c, a, b);
        nmod_poly_truncate(c, n);
        nmod_poly_sub_series(d, a, b, n);

        result = (nmod_poly_equal(d, c));
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

    TEST_FUNCTION_END(state);
}
