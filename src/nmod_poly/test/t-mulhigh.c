/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_mulhigh, state)
{
    int i, result;

    /* Compare with left truncated product of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        slong j, n;

        ulong m = n_randtest_not_zero(state);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        n = n_randint(state, 50);
        nmod_poly_randtest(b, state, n);
        nmod_poly_randtest(c, state, n);

        nmod_poly_mulhigh(a, b, c, n);
        nmod_poly_mul(b, b, c);
        for (j = 0; j < n; j++)
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

    /* Independent lengths and start, longer inputs, aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d;
        slong j, len1, len2, n;

        ulong m = n_randtest_not_zero(state);

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);

        if (n_randint(state, 10) == 0)
        {
            len1 = n_randint(state, 2000);
            len2 = n_randint(state, 2000);
        }
        else
        {
            len1 = n_randint(state, 100);
            len2 = n_randint(state, 100);
        }

        nmod_poly_randtest(b, state, len1);
        nmod_poly_randtest(c, state, len2);
        n = n_randint(state, len1 + len2 + 1);

        nmod_poly_mul(d, b, c);

        switch (n_randint(state, 3))
        {
            case 0:
                nmod_poly_mulhigh(a, b, c, n);
                break;
            case 1:
                nmod_poly_set(a, b);
                nmod_poly_mulhigh(a, a, c, n);
                break;
            default:
                nmod_poly_set(a, c);
                nmod_poly_mulhigh(a, b, a, n);
                break;
        }

        for (j = 0; j < n; j++)
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
            flint_printf("FAIL (lengths %wd, %wd, n = %wd):\n", len1, len2, n);
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
