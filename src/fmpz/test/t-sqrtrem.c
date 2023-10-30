/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_sqrtrem, state)
{
    int i, result;

    /* Comparison with mpz routines */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f, r, g;
        mpz_t mf, mf2, mr, mg;
        int aliasing;

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

        aliasing = n_randint(state, 3);

        if (aliasing == 0)
        {
            fmpz_sqrtrem(f, r, g);
        }
        else if (aliasing == 1)
        {
            fmpz_set(f, g);
            fmpz_sqrtrem(f, r, f);
        }
        else
        {
            fmpz_set(r, g);
            fmpz_sqrtrem(f, r, r);
        }

        mpz_sqrtrem(mf, mr, mg);

        fmpz_get_mpz(mf2, f);

        result = (mpz_cmp(mf2, mf) == 0) && _fmpz_is_canonical(f) && _fmpz_is_canonical(r);
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

    TEST_FUNCTION_END(state);
}
