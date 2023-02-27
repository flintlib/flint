/*
    Copyright (C) 2011 Sebastian Pancratz

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

    flint_printf("is_even/odd....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t g;

        fmpz_init(f);
        mpz_init(g);

        fmpz_randtest(f, state, 100);
        fmpz_get_mpz(g, f);

        result = (fmpz_is_even(f) == mpz_even_p(g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            gmp_printf("g = %Zd\n", g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        mpz_clear(g);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t g;

        fmpz_init(f);
        mpz_init(g);

        fmpz_randtest(f, state, 100);
        fmpz_get_mpz(g, f);

        result = (fmpz_is_odd(f) == mpz_odd_p(g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            gmp_printf("g = %Zd\n", g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        mpz_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
