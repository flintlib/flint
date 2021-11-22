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

    flint_printf("fdiv_ui....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t b;
        ulong x, r1, r2;

        fmpz_init(a);
        mpz_init(b);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(b, a);
        x = n_randtest_not_zero(state);

        r1 = fmpz_fdiv_ui(a, x);
        r2 = flint_mpz_fdiv_ui(b, x);

        result = (r1 == r2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("b = %Zd, x = %wu, r1 = %wu, r2 = %wu\n", b, x, r1, r2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        mpz_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
