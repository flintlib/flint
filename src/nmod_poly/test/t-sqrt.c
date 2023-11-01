/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_sqrt, state)
{
    int i;

    /* Test aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        int square1, square2;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);

        nmod_poly_randtest(a, state, 1 + n_randint(state, 50));

        if (n_randint(state, 2))
            nmod_poly_mul(a, a, a);

        square1 = nmod_poly_sqrt(b, a);
        square2 = nmod_poly_sqrt(a, a);

        if ((square1 != square2) || (square1 && !nmod_poly_equal(a, b)))
        {
            flint_printf("FAIL: aliasing:\n");
            flint_printf("square1 = %d, square2 = %d\n\n", square1, square2);
            flint_printf("a: "); nmod_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); nmod_poly_print(b); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* Test random squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        int square;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);
        nmod_poly_init(c, mod);

        nmod_poly_randtest(a, state, 1 + n_randint(state, 50));
        nmod_poly_mul(b, a, a);
        square = nmod_poly_sqrt(c, b);

        if (!square)
        {
            flint_printf("FAIL: square reported nonsquare:\n");
            flint_printf("a: "); nmod_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); nmod_poly_print(b); flint_printf("\n\n");
            flint_printf("c: "); nmod_poly_print(c); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mul(c, c, c);
        if (!nmod_poly_equal(c, b))
        {
            flint_printf("FAIL: sqrt(b)^2 != b:\n");
            flint_printf("a: "); nmod_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); nmod_poly_print(b); flint_printf("\n\n");
            flint_printf("c: "); nmod_poly_print(c); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Test "almost" squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        slong j;
        int square;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);
        nmod_poly_init(c, mod);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 50));
        nmod_poly_mul(b, a, a);

        j = n_randint(state, nmod_poly_length(b));
        b->coeffs[j] = n_randint(state, mod);
        _nmod_poly_normalise(b);

        square = nmod_poly_sqrt(c, b);

        if (square)
        {
            nmod_poly_mul(c, c, c);
            if (!nmod_poly_equal(c, b))
            {
                flint_printf("FAIL: sqrt(b)^2 != b:\n");
                flint_printf("a: "); nmod_poly_print(a); flint_printf("\n\n");
                flint_printf("b: "); nmod_poly_print(b); flint_printf("\n\n");
                flint_printf("c: "); nmod_poly_print(c); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
