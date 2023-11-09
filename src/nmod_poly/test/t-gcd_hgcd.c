/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_gcd_hgcd, state)
{
    int i, result;

    /*
       Find coprime polys, multiply by another poly
       and check the GCD is that poly
    */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, g;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(g, n);

        do {
            nmod_poly_randtest(a, state, n_randint(state, 1000));
            nmod_poly_randtest(b, state, n_randint(state, 1000));
            nmod_poly_gcd_hgcd(g, a, b);
        } while (g->length != 1);

        do {
            nmod_poly_randtest(c, state, n_randint(state, 1000));
        } while (c->length < 2);
        nmod_poly_make_monic(c, c);

        nmod_poly_mul(a, a, c);
        nmod_poly_mul(b, b, c);

        nmod_poly_gcd_hgcd(g, a, b);

        result = (nmod_poly_equal(g, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(g);
    }

    /* Check aliasing of a and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_randtest(a, state, n_randint(state, 1000));
        nmod_poly_randtest(b, state, n_randint(state, 1000));

        nmod_poly_gcd_hgcd(g, a, b);
        nmod_poly_gcd_hgcd(a, a, b);

        result = (nmod_poly_equal(a, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
    }

    /* Check aliasing of b and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_randtest(a, state, n_randint(state, 1000));
        nmod_poly_randtest(b, state, n_randint(state, 1000));

        nmod_poly_gcd_hgcd(g, a, b);
        nmod_poly_gcd_hgcd(b, a, b);

        result = (nmod_poly_equal(b, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
