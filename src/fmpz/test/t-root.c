/*
    Copyright (C) 2009, 2011 William Hart
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
    int i, result, exact;
    FLINT_TEST_INIT(state);

    flint_printf("root....");
    fflush(stdout);

    /* Comparison with mpz routines */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g;
        mpz_t mf, mf2, mg;
        slong n;

        fmpz_init(f);
        fmpz_init(g);

        mpz_init(mf);
        mpz_init(mf2);
        mpz_init(mg);

        n = n_randint(state, 20) + 1;
        
        fmpz_randtest(g, state, 200);
        if ((n & 1) == 0)
            fmpz_abs(g, g);
        fmpz_get_mpz(mg, g);

        fmpz_root(f, g, n);
        mpz_root(mf, mg, n);

        fmpz_get_mpz(mf2, f);

        result = (mpz_cmp(mf2, mf) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("mf = %Zd, mf2 = %Zd, mg = %Zd, root = %Md\n", mf, mf2, mg, n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        fmpz_clear(g);

        mpz_clear(mf);
        mpz_clear(mf2);
        mpz_clear(mg);
    }

    /* Exact powers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g, pow;
        slong n;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(pow);

        n = n_randint(state, 20) + 1;
        
        fmpz_randtest(g, state, 200);
        if ((n & 1) == 0)
            fmpz_abs(g, g);

        fmpz_pow_ui(pow, g, n);

        exact = fmpz_root(f, pow, n);

        result = (exact && fmpz_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("g = "); fmpz_print(g); flint_printf("\n");
            flint_printf("f = "); fmpz_print(f); flint_printf("\n");
            flint_printf("exact = %d, n = %wu\n", exact, n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(pow);
    }

    /* Not exact powers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, g, pow;
        slong n;

        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(pow);

        n = n_randint(state, 20) + 2;
        
        while (fmpz_cmp_ui(g, 2) < 0)
            fmpz_randtest_unsigned(g, state, 200);

        while (fmpz_is_zero(f))
            fmpz_randm(f, state, g);

        fmpz_pow_ui(pow, g, n);
        fmpz_add(pow, pow, f);

        if ((n & 1) != 0 && n_randint(state, 2) == 0)
        {
            fmpz_neg(g, g);
            fmpz_neg(pow, pow);
        }

        exact = fmpz_root(f, pow, n);

        result = (!exact && fmpz_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("g = "); fmpz_print(g); flint_printf("\n");
            flint_printf("f = "); fmpz_print(f); flint_printf("\n");
            flint_printf("exact = %d, n = %wu\n", exact, n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(pow);
    }

    /* Check aliasing of f and g */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t mf, mf2;
        slong n;

        fmpz_init(f);

        mpz_init(mf);
        mpz_init(mf2);

        n = n_randint(state, 20) + 1;
        
        fmpz_randtest(f, state, 200);
        if ((n & 1) == 0)
            fmpz_abs(f, f);
        fmpz_get_mpz(mf, f);

        fmpz_root(f, f, n);
        mpz_root(mf, mf, n);

        fmpz_get_mpz(mf2, f);

        result = (mpz_cmp(mf, mf2) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("mf = %Zd, mf2 = %Zd, root = %Md\n", mf, mf2, n);
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
