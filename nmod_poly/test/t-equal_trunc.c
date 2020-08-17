/*
    Copyright (C) 2009 William Hart

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
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("equal_trunc....");
    fflush(stdout);

    /* equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        ulong c, p;
        slong n, j;

        p = n_randtest(state);
	if (p == 0)
	   p++;

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        n = n_randint(state, 100);

        for (j = 0; j < n; j++)
        {
           c = nmod_poly_get_coeff_ui(a, j);
           nmod_poly_set_coeff_ui(b, j, c);
        }

        result = (nmod_poly_equal_trunc(a, b, n));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    /* unequal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        ulong c, p;
        slong m, n, j;

        p = n_randtest(state);
	if (p < 2)
	   p += 2;

        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        n = n_randint(state, 100) + 1;
        m = n_randint(state, n);

        for (j = 0; j < n; j++)
        {
           c = nmod_poly_get_coeff_ui(a, j);
           nmod_poly_set_coeff_ui(b, j, c);
        }
        c = nmod_poly_get_coeff_ui(b, m);
        c = n_addmod(c, 1, p);
        nmod_poly_set_coeff_ui(b, m, c);

        result = (!nmod_poly_equal_trunc(a, b, n));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            flint_printf("m = %wd\n", m);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
