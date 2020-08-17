/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);
    

    flint_printf("div_basecase....");
    fflush(stdout);

    /* Check result of div against divrem */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q, r, q2;

        mp_limb_t n;
        do
        {
            n = n_randtest_not_zero(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        nmod_poly_init(q2, n);
        
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length == 0);

        nmod_poly_divrem_basecase(q, r, a, b);
        nmod_poly_div_basecase(q2, a, b);

        result = (nmod_poly_equal(q, q2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(q2), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            abort();
        }
        
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
        nmod_poly_clear(q2);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q;

        mp_limb_t n;
        do
        {
            n = n_randtest(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_div_basecase(q, a, b);
        nmod_poly_div_basecase(a, a, b);

        result = (nmod_poly_equal(a, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q;

        mp_limb_t n;
        do
        {
            n = n_randtest(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_div_basecase(q, a, b);
        nmod_poly_div_basecase(b, a, b);

        result = (nmod_poly_equal(b, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
