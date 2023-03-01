/*
    Copyright (C) 2021 Daniel Schultz

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

    flint_printf("get_nmod....");
    fflush(stdout);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        ulong x, r1, r2;
        nmod_t mod;

        fmpz_init(a);

        fmpz_randtest(a, state, 1000);

        x = n_randtest_not_zero(state);

        nmod_init(&mod, x);

        r1 = fmpz_fdiv_ui(a, x);
        r2 = fmpz_get_nmod(a, mod);

        result = (r1 == r2);
        if (!result)
        {
            flint_printf("FAIL: i = %wd\n", i);
            flint_printf("a: "); fmpz_print(a); flint_printf("\n");
            flint_printf("x : %wu\n", x);
            flint_printf("r1: %wu\n", r1);
            flint_printf("r2: %wu\n", r2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
