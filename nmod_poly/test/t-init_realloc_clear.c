/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("init/init2/realloc/clear....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init2(a, n, n_randint(state, 100));
        nmod_poly_clear(a);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init2(a, n, n_randint(state, 100));
        nmod_poly_realloc(a, n_randint(state, 100));
        nmod_poly_realloc(a, n_randint(state, 100));
        nmod_poly_clear(a);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
