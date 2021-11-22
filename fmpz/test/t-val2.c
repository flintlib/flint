/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("val2....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t x;
        slong v1, v2;

        fmpz_init(x);

        /* Test special case */
        if (n_randint(state, 1000) == 0)
        {
            fmpz_zero(x);
            v1 = 0;
        }
        else
        {
            do {
                fmpz_randtest(x, state, 1000);
            } while (fmpz_is_even(x));

            v1 = n_randint(state, 1000);
            fmpz_mul_2exp(x, x, v1);
        }

        v2 = fmpz_val2(x);

        result = ((v1 == v2) == 1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("v1 = %wd  v2 = %wd\n", v1, v2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
