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
#include <limits.h>
#include "flint.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("divisible_si....");
    fflush(stdout);

    

    /* Compare with MPIR:  random */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong a;
        fmpz_t b;
        mpz_t d;
        int e, f;

        fmpz_init(b);
        mpz_init(d);

        a = z_randtest(state);
        if (a == WORD_MIN)
            a = 1;
        a = FLINT_ABS(a) + 1;
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(d, b);

        e = fmpz_divisible_si(b, a);
        f = flint_mpz_divisible_ui_p(d, a);

        result = (e == f);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = %wd, b = ", a), fmpz_print(b), flint_printf("\n");
            abort();
        }

        fmpz_clear(b);
        mpz_clear(d);
    }

    /* Compare with MPIR:  b a multiple of a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong a;
        fmpz_t b;
        mpz_t d;
        int e, f;

        fmpz_init(b);
        mpz_init(d);

        a = z_randtest(state);
        if (a == WORD_MIN)
            a = 1;
        a = FLINT_ABS(a) + 1;
        fmpz_randtest(b, state, 200);
        fmpz_mul_ui(b, b, a);

        fmpz_get_mpz(d, b);

        e = fmpz_divisible_si(b, a);
        f = flint_mpz_divisible_ui_p(d, a);

        result = (e == f && e == 1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = %wd, b = ", a), fmpz_print(b), flint_printf("\n");
            abort();
        }

        fmpz_clear(b);
        mpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

