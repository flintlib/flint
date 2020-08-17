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
    
    
    flint_printf("compose_horner....");
    fflush(stdout);

    /* Check (f(x-1))(x+1) == f */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, r, xp1, xm1;
        mp_limb_t n = n_randtest_not_zero(state);
        
        nmod_poly_init(a, n);
        nmod_poly_init(r, n);
        nmod_poly_init(xm1, n);
        nmod_poly_init(xp1, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        
        nmod_poly_set_coeff_ui(xm1, 1, 1);
        nmod_poly_set_coeff_ui(xm1, 0, n - 1);
        nmod_poly_set_coeff_ui(xp1, 1, 1);
        nmod_poly_set_coeff_ui(xp1, 0, 1);

        nmod_poly_compose_horner(r, a, xm1);
        nmod_poly_compose_horner(r, r, xp1);
        
        result = nmod_poly_equal(a, r);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(r);
        nmod_poly_clear(xm1);
        nmod_poly_clear(xp1);
    }

    /* Check a(c) + b(c) = (a + b)(c) */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, r1, r2;
        mp_limb_t n = n_randtest_not_zero(state);
        
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(r1, n);
        nmod_poly_init(r2, n);
        nmod_poly_randtest(a, state, n_randint(state, 30));
        nmod_poly_randtest(b, state, n_randint(state, 30));
        nmod_poly_randtest(c, state, n_randint(state, 10));
        
        nmod_poly_compose_horner(r1, a, c);
        nmod_poly_compose_horner(r2, b, c);
        nmod_poly_add(r1, r1, r2);

        nmod_poly_add(a, a, b);
        nmod_poly_compose_horner(r2, a, c);
        
        result = nmod_poly_equal(r1, r2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(r1), flint_printf("\n\n");
            nmod_poly_print(r2), flint_printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(r1);
        nmod_poly_clear(r2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
