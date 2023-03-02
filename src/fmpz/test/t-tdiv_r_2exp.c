/*
    Copyright (C) 2020 Daniel Schultz

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

    flint_printf("tdiv_r_2exp....");
    fflush(stdout);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f;
        ulong x;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randint(state, 200);

        fmpz_tdiv_r_2exp(b, a, x);
        mpz_tdiv_r_2exp(e, d, x);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0);
        if (!result)
        {
            flint_printf("FAIL 1:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, exp = %Mu\n", d, e, f, x);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t d, e, f;
        ulong x;

        fmpz_init(a);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randint(state, 200);

        fmpz_tdiv_r_2exp(a, a, x);
        mpz_tdiv_r_2exp(e, d, x);

        fmpz_get_mpz(f, a);

        result = (mpz_cmp(e, f) == 0);
        if (!result)
        {
            flint_printf("FAIL 2:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, exp = %Mu\n", d, e, f, x);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

