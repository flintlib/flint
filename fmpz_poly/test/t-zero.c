/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
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

    flint_printf("zero....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a;

        fmpz_poly_init(a);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 100);
        fmpz_poly_zero(a);

        result = (fmpz_poly_is_zero(a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
