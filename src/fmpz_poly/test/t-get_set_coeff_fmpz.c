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
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("get/set_coeff_fmpz....");
    fflush(stdout);

    

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a;
        fmpz_t x1, x2;
        slong coeff, len;

        fmpz_poly_init(a);
        fmpz_init(x1);
        fmpz_init(x2);
        len = n_randint(state, 100) + 1;

        for (j = 0; j < 1000; j++)
        {
            fmpz_randtest(x1, state, 200);
            coeff = n_randint(state, len);
            fmpz_poly_set_coeff_fmpz(a, coeff, x1);
            fmpz_poly_get_coeff_fmpz(x2, a, coeff);

            result = (fmpz_equal(x1, x2));
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("x1 = "), fmpz_print(x1), flint_printf("\n");
                flint_printf("x2 = "), fmpz_print(x2), flint_printf("\n");
                flint_printf("coeff = %wd, length = %wd\n", coeff, len);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(x1);
        fmpz_clear(x2);
        fmpz_poly_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
