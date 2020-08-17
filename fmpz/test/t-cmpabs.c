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

    flint_printf("cmpabs....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;

        fmpz_init(a);
        fmpz_init(b);

        fmpz_randtest(a, state, 200);

        fmpz_set(b, a);
        result = (fmpz_cmpabs(a, b));

        if (result != 0)
        {
            flint_printf("FAIL\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;
        int r1, r2;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(c);
        mpz_init(d);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(c, a);
        fmpz_get_mpz(d, b);

        r1 = fmpz_cmpabs(a, b);
        r2 = mpz_cmpabs(c, d);
        result = ((r1 == 0 && r2 == 0) || (r1 > 0 && r2 > 0)
                  || (r1 < 0 && r2 < 0));

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(c);
        mpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
