/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

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

    flint_printf("mul2_uiui....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f;
        ulong x, y;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randtest(state);
        y = n_randtest(state);

        fmpz_mul2_uiui(b, a, x, y);
        flint_mpz_mul_ui(e, d, x);
        flint_mpz_mul_ui(e, e, y);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, x = %Mu, y = %Mu\n",
                d, e, f, x, y);
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
        ulong x, y;

        fmpz_init(a);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randtest(state);
        y = n_randtest(state);

        fmpz_mul2_uiui(a, a, x, y);
        flint_mpz_mul_ui(e, d, x);
        flint_mpz_mul_ui(e, e, y);

        fmpz_get_mpz(f, a);

        result = (mpz_cmp(e, f) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, x = %Mu, y = %Mu\n",
                d, e, f, x, y);
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
    return 0;
}
