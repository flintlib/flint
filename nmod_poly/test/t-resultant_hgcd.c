/*
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
    

    flint_printf("resultant_hgcd....");
    fflush(stdout);

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g;
        mp_limb_t x, y;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(f, n);
        nmod_poly_init(g, n);
        
        nmod_poly_randtest(f, state, n_randint(state, 200));
        nmod_poly_randtest(g, state, n_randint(state, 200));

        x = nmod_poly_resultant_hgcd(f, g);
        y = nmod_poly_resultant_hgcd(g, f);

        if ((nmod_poly_degree(f) * nmod_poly_degree(g)) % 2)
            y = nmod_neg(y, f->mod);

        result = (x == y);
        if (!result)
        {
            flint_printf("FAIL (res(f, g) == (-1)^(deg f deg g) res(g, f)):\n");
            nmod_poly_print(f), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("x = %wu\n", x);
            flint_printf("y = %wu\n", y);
            flint_printf("n = %wu\n", n);
            abort();
        }
        
        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h;
        mp_limb_t x, y, z;
        mp_limb_t n;

        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(f, n);
        nmod_poly_init(g, n);
        nmod_poly_init(h, n);
        
        nmod_poly_randtest(f, state, n_randint(state, 600));
        nmod_poly_randtest(g, state, n_randint(state, 600));
        nmod_poly_randtest(h, state, n_randint(state, 600));

        y = nmod_poly_resultant_hgcd(f, g);
        z = nmod_poly_resultant_hgcd(h, g);
        y = nmod_mul(y, z, f->mod);
        nmod_poly_mul(f, f, h);
        x = nmod_poly_resultant_hgcd(f, g);

        result = (x == y);
        if (!result)
        {
            flint_printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            nmod_poly_print(f), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            nmod_poly_print(h), flint_printf("\n\n");
            flint_printf("x = %wu\n", x);
            flint_printf("y = %wu\n", y);
            flint_printf("z = %wd\n", z);
            flint_printf("n = %wu\n", n);
            abort();
        }
        
        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
