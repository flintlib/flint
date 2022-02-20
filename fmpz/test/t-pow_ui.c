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

    flint_printf("pow_ui....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g, tst;
        mpz_t mf, mg;
        ulong exp;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(tst);

        mpz_init(mf);
        mpz_init(mg);

        fmpz_randtest(g, state, 200);

        fmpz_get_mpz(mg, g);
        exp = n_randint(state, 70);

        fmpz_pow_ui(f, g, exp);
        flint_mpz_pow_ui(mf, mg, exp);

        fmpz_set_mpz(tst, mf);

        result = fmpz_equal(f, tst);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "); fmpz_print(f); flint_printf(", ");
            flint_printf("g = "); fmpz_print(g); flint_printf(", ");
            flint_printf("exp = %wu\n", exp);
            flint_printf("Correct output via GMP: "); fmpz_print(tst); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(tst);

        mpz_clear(mf);
        mpz_clear(mg);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g, tst;
        mpz_t mf;
        ulong exp;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(tst);

        mpz_init(mf);

        fmpz_randtest(g, state, 200);
        fmpz_set(f, g);

        fmpz_get_mpz(mf, f);
        exp = n_randint(state, 70);

        fmpz_pow_ui(f, f, exp);
        flint_mpz_pow_ui(mf, mf, exp);

        fmpz_set_mpz(tst, mf);

        result = fmpz_equal(f, tst);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("(aliasing)\n");
            flint_printf("f = "); fmpz_print(f); flint_printf(", ");
            flint_printf("g = "); fmpz_print(g); flint_printf(", ");
            flint_printf("exp = %wu\n", exp);
            flint_printf("Correct output via GMP: "); fmpz_print(tst); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(tst);

        mpz_clear(mf);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
