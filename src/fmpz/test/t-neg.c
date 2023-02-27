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

    flint_printf("neg....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(c);
        mpz_init(d);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(c, a);

        fmpz_neg(b, a);
        mpz_neg(c, c);

        fmpz_get_mpz(d, b);

        result = (mpz_cmp(c, d) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(c);
        mpz_clear(d);
    }

    /* Check aliasing */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t c, d;

        fmpz_init(a);

        mpz_init(c);
        mpz_init(d);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(c, a);

        fmpz_neg(a, a);
        mpz_neg(c, c);

        fmpz_get_mpz(d, a);

        result = (mpz_cmp(c, d) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);

        mpz_clear(c);
        mpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
