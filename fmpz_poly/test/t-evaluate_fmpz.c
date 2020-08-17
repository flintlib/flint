/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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
#include "fmpz_poly.h"
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
        fmpz_t a, b;
        fmpz_poly_t f;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpz_randtest(a, state, 100);

        fmpz_poly_evaluate_fmpz(b, f, a);
        fmpz_poly_evaluate_fmpz(a, f, a);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_poly_clear(f);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
