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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;

    FLINT_TEST_INIT(state);
    

    flint_printf("set_uiui....");
    fflush(stdout);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        ulong hi, lo;

        hi = n_randtest(state);
        lo = n_randtest(state);

        fmpz_init(a);
        fmpz_init(b);

        fmpz_set_ui(a, hi);
        fmpz_mul_2exp(a, a, FLINT_BITS);
        fmpz_add_ui(a, a, lo);

        fmpz_set_uiui(b, hi, lo);

        result = fmpz_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("hi = %wu\n", hi);
            flint_printf("lo = %wu\n", lo);
            flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
