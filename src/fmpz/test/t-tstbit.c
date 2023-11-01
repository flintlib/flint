/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_tstbit, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        int k, l;
        ulong j;
        fmpz_t a;
        mpz_t b;

        fmpz_init(a);
        mpz_init(b);

        fmpz_randtest(a, state, 2 * FLINT_BITS);
        fmpz_get_mpz(b, a);
        j = n_randint(state, 3 * FLINT_BITS);

        k = fmpz_tstbit(a, j);
        l = mpz_tstbit(b, j);

        result = (k == l);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd, j = %Mu k = %d, l = %d\n", b, j, k, l);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        mpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
