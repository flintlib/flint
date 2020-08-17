/*
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

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("get_coeff_ptr....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t A;
        fmpz_t a;
        slong n = n_randint(state, 100);

        fmpz_poly_init(A);
        fmpz_poly_randtest(A, state, n_randint(state, 100), 100);
        fmpz_init(a);

        fmpz_poly_get_coeff_fmpz(a, A, n);

        result = n < fmpz_poly_length(A) ? 
                     fmpz_equal(a, fmpz_poly_get_coeff_ptr(A, n)) : 
                     fmpz_poly_get_coeff_ptr(A, n) == NULL;
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A = "), fmpz_poly_print(A), flint_printf("\n\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n\n");
            flint_printf("n = %wd\n\n", n);
            abort();
        }

        fmpz_poly_clear(A);
        fmpz_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
