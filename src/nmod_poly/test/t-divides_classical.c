/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_divides_classical, state)
{
    int i, result;

    /* Random polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q, prod;
        int divides;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_init(prod, n);

        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        divides = nmod_poly_divides_classical(q, a, b);
        nmod_poly_mul(prod, q, b);

        result = ((divides && nmod_poly_equal(a, prod)) ||
		 (!divides && nmod_poly_is_zero(q) && !nmod_poly_equal(prod, a)));
        if (!result)
        {
            flint_printf("FAIL:\n");
	    flint_printf("divides = %d\n", divides);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
	    nmod_poly_print(prod), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(prod);
    }

    /* Random divisible polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q, prod;
        int divides;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_init(prod, n);

        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        nmod_poly_mul(a, b, a);

        divides = nmod_poly_divides_classical(q, a, b);
        nmod_poly_mul(prod, q, b);

        result = (divides && nmod_poly_equal(a, prod));
        if (!result)
        {
            flint_printf("FAIL:\n");
	    flint_printf("divides = %d\n", divides);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
	    nmod_poly_print(prod), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(prod);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q;
        int divides1, divides2;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);

        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        divides1 = nmod_poly_divides_classical(q, a, b);
        divides2 = nmod_poly_divides_classical(a, a, b);

        result = (divides1 == divides2 && nmod_poly_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
	    flint_printf("divides1 = %d, divides2 = %d\n", divides1, divides2);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q;
        int divides1, divides2;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);

        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));

        divides1 = nmod_poly_divides_classical(q, a, b);
        divides2 = nmod_poly_divides_classical(b, a, b);

        result = (divides1 == divides2 && nmod_poly_equal(q, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
	    flint_printf("divides1 = %d, divides2 = %d\n", divides1, divides2);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}
