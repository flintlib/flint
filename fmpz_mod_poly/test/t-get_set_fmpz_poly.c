/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
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

    flint_printf("get/set_fmpz_poly....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b;
        fmpz_poly_t c;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_poly_init(c);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));

        fmpz_mod_poly_get_fmpz_poly(c, a);
        fmpz_mod_poly_set_fmpz_poly(b, c);

        result = fmpz_mod_poly_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_mod_poly_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n");
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_poly_clear(c);

        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

