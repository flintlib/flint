/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("discriminant....");
    fflush(stdout);

    /* Check disc(fg) == disc(f) * disc(g) * res(f, g)^2 */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t f, g, h;
        fmpz_t x, y, z, r;
        fmpz_t n;

        fmpz_init(n);
        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);
        fmpz_init(r);

        fmpz_set_ui(n, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(f, n);
        fmpz_mod_poly_init(g, n);
        fmpz_mod_poly_init(h, n);
        
        do {
           fmpz_mod_poly_randtest(f, state, n_randint(state, 200));
        } while (f->length < 2);
        do {
           fmpz_mod_poly_randtest(g, state, n_randint(state, 200));
        } while (g->length < 2);

        fmpz_mod_poly_discriminant(y, f);
        fmpz_mod_poly_discriminant(z, g);
        fmpz_mul(y, y, z);
        fmpz_mod(y, y, n);
        fmpz_mod_poly_resultant(r, f, g);
        fmpz_mul(r, r, r);
        fmpz_mod(r, r, n);
        fmpz_mul(y, y, r);
        fmpz_mod(y, y, n);
        fmpz_mod_poly_mul(h, f, g);
        fmpz_mod_poly_discriminant(x, h);

        result = (fmpz_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL (disc(fg) == res(f, g)^2 * disc(f) * disc(g):\n");
            fmpz_mod_poly_print(f), flint_printf("\n\n");
            fmpz_mod_poly_print(g), flint_printf("\n\n");
            fmpz_mod_poly_print(h), flint_printf("\n\n");
            flint_printf("x = "); fmpz_print(x); printf("\n");
            flint_printf("y = "); fmpz_print(y); printf("\n");
            flint_printf("z = "); fmpz_print(z); printf("\n");
            flint_printf("n = "); fmpz_print(n); printf("\n");
            abort();
        }
        
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(h);

        fmpz_clear(n);
        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
        fmpz_clear(r);
    }

    /* Check disc(f) == 0 for length < 2 */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t f;
        fmpz_t y;
        fmpz_t n;

        fmpz_init(y);
        fmpz_init(n);

        fmpz_set_ui(n, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(f, n);
        
        fmpz_mod_poly_randtest(f, state, 1);
        
        fmpz_mod_poly_discriminant(y, f);
        
        result = fmpz_is_zero(y);
        if (!result)
        {
            flint_printf("FAIL disc(f) == 0 for len f < 2:\n");
            fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("y = "); fmpz_print(y); printf("\n");
            abort();
        }
        
        fmpz_clear(n);
        fmpz_clear(y);

        fmpz_mod_poly_clear(f);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
