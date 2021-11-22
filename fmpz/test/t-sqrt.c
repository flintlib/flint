/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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

    flint_printf("sqrt....");
    fflush(stdout);

    

    /* Comparison with mpz routines */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g;
        mpz_t mf, mf2, mg;

        fmpz_init(f);
        fmpz_init(g);

        mpz_init(mf);
        mpz_init(mf2);
        mpz_init(mg);

        fmpz_randtest(g, state, 200);
        fmpz_abs(g, g);
        fmpz_get_mpz(mg, g);

        fmpz_sqrt(f, g);
        mpz_sqrt(mf, mg);

        fmpz_get_mpz(mf2, f);

        result = (mpz_cmp(mf2, mf) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("mf = %Zd, mf2 = %Zd, mg = %Zd\n", mf, mf2, mg);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        fmpz_clear(g);

        mpz_clear(mf);
        mpz_clear(mf2);
        mpz_clear(mg);
    }

    /* Check aliasing of f and g */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t mf, mf2;

        fmpz_init(f);

        mpz_init(mf);
        mpz_init(mf2);

        fmpz_randtest(f, state, 200);
        fmpz_abs(f, f);
        fmpz_get_mpz(mf, f);

        fmpz_sqrt(f, f);
        mpz_sqrt(mf, mf);

        fmpz_get_mpz(mf2, f);

        result = (mpz_cmp(mf, mf2) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("mf = %Zd, mf2 = %Zd\n", mf, mf2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);

        mpz_clear(mf);
        mpz_clear(mf2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
