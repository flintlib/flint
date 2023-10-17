/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_kronecker, state)
{
    slong i, j;

    for (i = 0; i < 3000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, n;
        mpz_t aa, nn;

        fmpz_init(a);
        fmpz_init(n);
        mpz_init(aa);
        mpz_init(nn);

        for (j = 0; j < 100; j++)
        {
            fmpz_randtest(a, state, 150);
            fmpz_randtest(n, state, 150);

            fmpz_get_mpz(aa, a);
            fmpz_get_mpz(nn, n);

            if (mpz_kronecker(aa, nn) != fmpz_kronecker(a, n))
            {
                flint_printf("FAIL:\n");
                gmp_printf("a = %Zd, n = %Zd\n", aa, nn);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(a);
        fmpz_clear(n);
        mpz_clear(aa);
        mpz_clear(nn);
    }

    TEST_FUNCTION_END(state);
}
