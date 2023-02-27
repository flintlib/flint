/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Sebastian Pancratz

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

    flint_printf("tstbit....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        int k, l;
        ulong j;
        fmpz_t a;
        mpz_t b;

        fmpz_init(a);
        mpz_init(b);

        fmpz_randtest(a, state, 2 * FLINT_BITS);
        fmpz_get_mpz(b, a);
        j = n_randint(state, 3 * FLINT_BITS);

        k = fmpz_tstbit(a, j);
        l = mpz_tstbit(b, j);

        result = (k == l);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd, j = %Mu k = %d, l = %d\n", b, j, k, l);
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
