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

    flint_printf("sqrtrem....");
    fflush(stdout);

    

    /* Comparison with mpz routines */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, r, g;
        mpz_t mf, mf2, mr, mg;

        fmpz_init(f);
        fmpz_init(r);
        fmpz_init(g);

        mpz_init(mf);
        mpz_init(mf2);
        mpz_init(mr);
        mpz_init(mg);

        fmpz_randtest(g, state, 200);
        fmpz_abs(g, g);
        fmpz_get_mpz(mg, g);

        fmpz_sqrtrem(f, r, g);
        mpz_sqrtrem(mf, mr, mg);

        fmpz_get_mpz(mf2, f);

        result = (mpz_cmp(mf2, mf) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("mf = %Zd, mf2 = %Zd, mr = %Zd, mg = %Zd\n", mf, mf2,
                       mr, mg);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        fmpz_clear(r);
        fmpz_clear(g);

        mpz_clear(mf);
        mpz_clear(mf2);
        mpz_clear(mr);
        mpz_clear(mg);
    }

    /* Check aliasing of r and g */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, f;
        mpz_t ma, mf, mf2;

        fmpz_init(a);
        fmpz_init(f);

        mpz_init(ma);
        mpz_init(mf);
        mpz_init(mf2);

        fmpz_randtest(a, state, 200);
        fmpz_abs(a, a);
        fmpz_get_mpz(ma, a);

        fmpz_sqrtrem(f, a, a);
        mpz_sqrtrem(mf, ma, ma);

        fmpz_get_mpz(mf2, f);

        result = (mpz_cmp(mf, mf2) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("ma = %Zd, mf = %Zd, mf2 = %Zd\n", ma, mf, mf2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(f);

        mpz_clear(ma);
        mpz_clear(mf);
        mpz_clear(mf2);
    }

    /* Check aliasing of f and g */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t r, f;
        mpz_t mr, mf, mf2;

        fmpz_init(r);
        fmpz_init(f);

        mpz_init(mr);
        mpz_init(mf);
        mpz_init(mf2);

        fmpz_randtest(f, state, 200);
        fmpz_abs(f, f);
        fmpz_get_mpz(mf, f);

        fmpz_sqrtrem(f, r, f);
        mpz_sqrtrem(mf, mr, mf);

        fmpz_get_mpz(mf2, f);

        result = (mpz_cmp(mf, mf2) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("mr = %Zd, mf = %Zd, mf2 = %Zd\n", mr, mf, mf2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(r);
        fmpz_clear(f);

        mpz_clear(mr);
        mpz_clear(mf);
        mpz_clear(mf2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
