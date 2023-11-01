/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_div_series, state)
{
    int i, result;

    /* Check A/B * B = A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t q, a, b, prod;
        slong m;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(prod, n);
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);

        nmod_poly_randtest(a, state, n_randint(state, 2000));
        do nmod_poly_randtest(b, state, n_randint(state, 2000));
        while (b->length == 0 || b->coeffs[0] == 0);

        m = n_randint(state, 2000) + 1;

        nmod_poly_div_series(q, a, b, m);
        nmod_poly_mullow(prod, q, b, m);
        nmod_poly_truncate(a, m);

        result = (nmod_poly_equal(a, prod));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(prod), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(q);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(prod);
    }

    /* Check aliasing of q and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t q, a, b;
        slong m;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(q, n);
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);

        nmod_poly_randtest(a, state, n_randint(state, 1000));
        do nmod_poly_randtest(b, state, n_randint(state, 1000));
        while (b->length == 0 || b->coeffs[0] == 0);

        m = n_randint(state, 1000) + 1;

        nmod_poly_div_series(q, a, b, m);
        nmod_poly_div_series(a, a, b, m);

        result = (nmod_poly_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("n = %wd, m = %wd\n", n, m);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(q);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* Check aliasing of q and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t q, a, b;
        slong m;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(q, n);
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);

        nmod_poly_randtest(a, state, n_randint(state, 1000));
        do nmod_poly_randtest(b, state, n_randint(state, 1000));
        while (b->length == 0 || b->coeffs[0] == 0);

        m = n_randint(state, 1000) + 1;

        nmod_poly_div_series(q, a, b, m);
        nmod_poly_div_series(b, a, b, m);

        result = (nmod_poly_equal(q, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("n = %wd, m = %wd\n", n, m);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(q);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
