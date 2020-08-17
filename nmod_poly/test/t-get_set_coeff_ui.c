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
    int i, result;
    ulong j;
    FLINT_TEST_INIT(state);
    

    flint_printf("get/set_coeff_ui....");
    fflush(stdout);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);
        mp_limb_t c1 = n_randtest(state), c2;
        
        j = n_randint(state, 100);

        nmod_poly_init(a, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        nmod_poly_set_coeff_ui(a, j, c1);
        c2 = nmod_poly_get_coeff_ui(a, j);
        
        result = (c2 == c1 % n);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("j = %wu, c1 = %wu, c2 = %wu, n = %wu\n", j, c1, c2, a->mod.n);
            nmod_poly_print(a), flint_printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
