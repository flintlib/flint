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

    flint_printf("mul....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(c, state, n_randint(state, 50));

        fmpz_mod_poly_mul(a, b, c);
        fmpz_mod_poly_mul(b, b, c);

        result = (fmpz_mod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a), flint_printf("\n\n");
            fmpz_mod_poly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_clear(p);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 50));
        fmpz_mod_poly_randtest(c, state, n_randint(state, 50));

        fmpz_mod_poly_mul(a, b, c);
        fmpz_mod_poly_mul(c, b, c);

        result = (fmpz_mod_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a), flint_printf("\n\n");
            fmpz_mod_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_clear(p);
    }

    /* Check (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a1, a2, b, c, d;
        fmpz_t p;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a1, p);
        fmpz_mod_poly_init(a2, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(d, p);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(c, state, n_randint(state, 100));
        fmpz_mod_poly_randtest(d, state, n_randint(state, 100));

        fmpz_mod_poly_mul(a1, b, c);
        fmpz_mod_poly_mul(a2, b, d);
        fmpz_mod_poly_add(a1, a1, a2);

        fmpz_mod_poly_add(c, c, d);
        fmpz_mod_poly_mul(a2, b, c);

        result = (fmpz_mod_poly_equal(a1, a2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a1), flint_printf("\n\n");
            fmpz_mod_poly_print(a2), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(a1);
        fmpz_mod_poly_clear(a2);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(d);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
