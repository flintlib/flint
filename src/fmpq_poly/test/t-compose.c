/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_compose, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpq_poly_compose(f, g, h);
        fmpq_poly_compose(g, g, h);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (aliasing 1):\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpq_poly_compose(f, g, h);
        fmpq_poly_compose(h, g, h);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(h) ? 0 : 2;
        result = (fmpq_poly_equal(f, h) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (aliasing 2):\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            fmpq_poly_debug(h), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Compare with the naive method for g(h(t)) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h, s, t, u;
        fmpq_t c;
        slong k;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_init(u);
        fmpq_init(c);
        fmpq_poly_randtest(g, state, n_randint(state, 20), 65);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 65);

        fmpq_poly_zero(s);
        fmpq_poly_set_ui(t, UWORD(1));
        for (k = WORD(0); k < g->length; k++)
        {
            fmpq_poly_get_coeff_fmpq(c, g, k);
            fmpq_poly_scalar_mul_fmpq(u, t, c);
            fmpq_poly_add(s, s, u);
            fmpq_poly_mul(t, t, h);
        }

        fmpq_poly_compose(f, g, h);

        result = (fmpq_poly_equal(f, s));
        if (!result)
        {
            flint_printf("FAIL (compare with naive):\n");
            flint_printf("g = "), fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("h = "), fmpq_poly_debug(h), flint_printf("\n\n");
            flint_printf("f = "), fmpq_poly_debug(f), flint_printf("\n\n");
            flint_printf("s = "), fmpq_poly_debug(s), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
        fmpq_poly_clear(u);
        fmpq_clear(c);
    }

    TEST_FUNCTION_END(state);
}
