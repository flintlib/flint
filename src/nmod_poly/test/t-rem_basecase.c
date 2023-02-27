/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Sebastian Pancratz

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
    

    flint_printf("rem_basecase....");
    fflush(stdout);

    /* Check result of divrem */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q0, r0, r;

        mp_limb_t n;
        do
        {
            n = n_randtest_not_zero(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q0, n);
        nmod_poly_init(r0, n);
        nmod_poly_init(r, n);
        
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_divrem_basecase(q0, r0, a, b);
        nmod_poly_rem_basecase(r, a, b);

        result = (nmod_poly_equal(r0, r));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q0), flint_printf("\n\n");
            nmod_poly_print(r0), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }
        
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q0);
        nmod_poly_clear(r0);
        nmod_poly_clear(r);
    }

    /* Check aliasing of a and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, r;

        mp_limb_t n;
        do
        {
            n = n_randtest(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(r, n);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_rem_basecase(r, a, b);
        nmod_poly_rem_basecase(a, a, b);

        result = (nmod_poly_equal(a, r));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r);
    }

    /* Check aliasing of b and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, r;

        mp_limb_t n;
        do
        {
            n = n_randtest(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(r, n);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_rem_basecase(r, a, b);
        nmod_poly_rem_basecase(b, a, b);

        result = (nmod_poly_equal(b, r));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

