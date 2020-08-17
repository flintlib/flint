/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
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
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("evaluate_fmpz....");
    fflush(stdout);

    

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, p;
        fmpz_mod_poly_t f;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_randm(a, state, p);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 100));

        fmpz_mod_poly_evaluate_fmpz(b, f, a);
        fmpz_mod_poly_evaluate_fmpz(a, f, a);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(f), flint_printf("\n\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(f);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(p);
    }

    /* Check that the result agrees with Z[X] */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, p;
        fmpz_mod_poly_t f;
        fmpz_poly_t g;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_mod_poly_init(f, p);
        fmpz_poly_init(g);
        fmpz_randm(a, state, p);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 100));
        fmpz_mod_poly_get_fmpz_poly(g, f);

        fmpz_mod_poly_evaluate_fmpz(b, f, a);
        fmpz_poly_evaluate_fmpz(c, g, a);
        fmpz_mod(c, c, p);

        result = (fmpz_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (cmp with fmpz_poly):\n");
            fmpz_mod_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(g), flint_printf("\n\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            fmpz_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
        fmpz_mod_poly_clear(f);
        fmpz_poly_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

