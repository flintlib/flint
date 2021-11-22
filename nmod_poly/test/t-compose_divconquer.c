/*
    Copyright (C) 2010 William Hart

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
    int i, result = 1;
    FLINT_TEST_INIT(state);
    
    
    flint_printf("compose_divconquer....");
    fflush(stdout);

    /* Compare aliasing */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, r1;
        mp_limb_t n = n_randtest_not_zero(state);
        
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(r1, n);
        nmod_poly_randtest(a, state, n_randint(state, 30));
        nmod_poly_randtest(b, state, n_randint(state, 15));
        
        nmod_poly_compose_divconquer(r1, a, b);
        nmod_poly_compose_divconquer(a, a, b);
        
        result = nmod_poly_equal(r1, a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            nmod_poly_print(r1), flint_printf("\n\n");
            nmod_poly_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r1);
    }
    
    /* Compare other aliasing */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, r1;
        mp_limb_t n = n_randtest_not_zero(state);
        
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(r1, n);
        nmod_poly_randtest(a, state, n_randint(state, 30));
        nmod_poly_randtest(b, state, n_randint(state, 15));
        
        nmod_poly_compose_divconquer(r1, a, b);
        nmod_poly_compose_divconquer(b, a, b);
        
        result = nmod_poly_equal(r1, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            nmod_poly_print(r1), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r1);
    }
    
    /* Compare with compose_horner */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, r1, r2;
        mp_limb_t n = n_randtest_not_zero(state);
        
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(r1, n);
        nmod_poly_init(r2, n);
        nmod_poly_randtest(a, state, n_randint(state, 30));
        nmod_poly_randtest(b, state, n_randint(state, 15));
        
        nmod_poly_compose_divconquer(r1, a, b);
        nmod_poly_compose_horner(r2, a, b);
        
        result = nmod_poly_equal(r1, r2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            nmod_poly_print(r1), flint_printf("\n\n");
            nmod_poly_print(r2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r1);
        nmod_poly_clear(r2);
    }
    
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
