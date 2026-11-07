/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_mul_KS, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        nmod_poly_mul_KS(a, b, c);
        nmod_poly_mul_KS(b, b, c);

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
        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        nmod_poly_mul_KS(a, b, c);
        nmod_poly_mul_KS(c, b, c);

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

    /* Compare with mul_classical */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a1, a2, b, c;
        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a1, n);
        nmod_poly_init(a2, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        nmod_poly_mul_classical(a1, b, c);
        nmod_poly_mul_KS(a2, b, c);

        result = (nmod_poly_equal(a1, a2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a1), flint_printf("\n\n");
            nmod_poly_print(a2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a1);
        nmod_poly_clear(a2);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /*
        All coefficients n - 1, so that the coefficients of the product reach
        their largest values, with n of nu > FLINT_BITS / 2 bits (two- and
        three-word coefficients), sometimes 2^nu - 1; against
        nmod_poly_mul_classical
    */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a1, a2, b, c;
        flint_bitcnt_t nu = FLINT_BITS / 2 + 1 + n_randint(state, FLINT_BITS / 2);
        ulong n = n_randint(state, 4) ? n_randbits(state, nu)
                    : ((nu == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << nu) - 1);
        slong j, len1 = 1 + n_randint(state, 64), len2 = 1 + n_randint(state, 64);

        nmod_poly_init(a1, n);
        nmod_poly_init(a2, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        for (j = 0; j < len1; j++)
            nmod_poly_set_coeff_ui(b, j, n - 1);
        for (j = 0; j < len2; j++)
            nmod_poly_set_coeff_ui(c, j, n - 1);

        nmod_poly_mul_classical(a1, b, c);
        nmod_poly_mul_KS(a2, b, c);

        if (!nmod_poly_equal(a1, a2))
            TEST_FUNCTION_FAIL("all coefficients n - 1\nn = %wu, len1 = %wd, len2 = %wd\n"
                               "a1 = %{nmod_poly}\na2 = %{nmod_poly}\n", n, len1, len2, a1, a2);

        nmod_poly_clear(a1);
        nmod_poly_clear(a2);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
