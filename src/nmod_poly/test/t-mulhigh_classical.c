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

TEST_FUNCTION_START(nmod_poly_mulhigh_classical, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        slong j, start;

        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));
        start = n_randint(state, 50);

        nmod_poly_mulhigh_classical(a, b, c, start);
        nmod_poly_mulhigh_classical(b, b, c, start);

        for (j = 0; j < start; j++)
        {
            if (j < a->length)
                a->coeffs[j] = 0;
            if (j < b->length)
                b->coeffs[j] = 0;
        }
        _nmod_poly_normalise(a);
        _nmod_poly_normalise(b);

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        slong j, start;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));
        start = n_randint(state, 50);

        nmod_poly_mulhigh_classical(a, b, c, start);
        nmod_poly_mulhigh_classical(c, b, c, start);

        for (j = 0; j < start; j++)
        {
            if (j < a->length)
                a->coeffs[j] = 0;
            if (j < c->length)
                c->coeffs[j] = 0;
        }
        _nmod_poly_normalise(a);
        _nmod_poly_normalise(c);

        result = (nmod_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Compare with mul_basecase */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        slong j, start;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));
        start = n_randint(state, 50);

        nmod_poly_mul_classical(a, b, c);
        nmod_poly_mulhigh_classical(d, b, c, start);

        for (j = 0; j < start; j++)
        {
            if (j < a->length)
                a->coeffs[j] = 0;
            if (j < d->length)
                d->coeffs[j] = 0;
        }
        _nmod_poly_normalise(a);
        _nmod_poly_normalise(d);

        result = (nmod_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            nmod_poly_print(a), flint_printf("\n\n");
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
