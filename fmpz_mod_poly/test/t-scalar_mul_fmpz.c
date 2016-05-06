/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart

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

    flint_printf("scalar_mul_fmpz....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t n, p;
        fmpz_mod_poly_t a, b;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_init(n);
        fmpz_randtest(n, state, 200);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100));

        fmpz_mod_poly_scalar_mul_fmpz(b, a, n);
        fmpz_mod_poly_scalar_mul_fmpz(a, a, n);

        result = (fmpz_mod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a), flint_printf("\n\n");
            fmpz_mod_poly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(n);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

