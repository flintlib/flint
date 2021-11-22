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

    flint_printf("combit....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong j;
        fmpz_t a;
        mpz_t b, c;

        fmpz_init(a);
        mpz_init(b);
        mpz_init(c);

        fmpz_randtest(a, state, 2 * FLINT_BITS);
        fmpz_get_mpz(b, a);
        j = n_randint(state, 3 * FLINT_BITS);

        fmpz_combit(a, j);
        mpz_combit(b, j);
        fmpz_get_mpz(c, a);

        result = (mpz_cmp(b, c) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            gmp_printf("b = %Zd\n", b);
            flint_printf("j = %wd\n", j);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        mpz_clear(c);
        mpz_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
