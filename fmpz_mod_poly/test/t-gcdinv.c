/*
    Copyright (C) 2012 Sebastian Pancratz

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
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("gcdinv....");
    fflush(stdout);

    

    /* Generic case, most likely co-prime arguments ******************************/

    /* Check s*a == g mod b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, g, s, u;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(u, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        do 
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 2);

        fmpz_mod_poly_gcdinv(g, s, a, b);
        fmpz_mod_poly_mul(u, s, a);
        fmpz_mod_poly_rem(u, u, b);
        
        result = (fmpz_mod_poly_equal(g, u) || fmpz_mod_poly_is_zero(g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("g = "), fmpz_mod_poly_print(g), flint_printf("\n\n");
            flint_printf("s = "), fmpz_mod_poly_print(s), flint_printf("\n\n");
            flint_printf("u = "), fmpz_mod_poly_print(u), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(u);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Check s*a == g mod b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, g, s, u;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(g, p);
        fmpz_mod_poly_init(s, p);
        fmpz_mod_poly_init(u, p);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));
        do 
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        while (b->length < 2);
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_mul(a, f, a);
        fmpz_mod_poly_mul(b, f, b);

        fmpz_mod_poly_gcdinv(g, s, a, b);
        fmpz_mod_poly_mul(u, s, a);
        fmpz_mod_poly_rem(u, u, b);

        result = (fmpz_mod_poly_equal(u, g) || fmpz_mod_poly_is_zero(g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n\n");
            flint_printf("f = "), fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpz_mod_poly_print(g), flint_printf("\n\n");
            flint_printf("s = "), fmpz_mod_poly_print(s), flint_printf("\n\n");
            flint_printf("u = "), fmpz_mod_poly_print(u), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(s);
        fmpz_mod_poly_clear(u);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

