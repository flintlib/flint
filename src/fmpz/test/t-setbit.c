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

TEST_FUNCTION_START(fmpz_setbit, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong j;
        fmpz_t a, b, c;
        mpz_t z;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        mpz_init(z);

        fmpz_randtest(a, state, 2 * FLINT_BITS);
        fmpz_set(b, a);
        fmpz_get_mpz(z, b);
        j = n_randint(state, 3 * FLINT_BITS);

        fmpz_setbit(b, j);
        mpz_setbit(z, j);
        fmpz_set_mpz(c, z);

        result = (fmpz_equal(b, c)) && _fmpz_is_canonical(b);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            gmp_printf("z = %Zd\n", z);
            flint_printf("j = %wd\n", j);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        mpz_clear(z);
    }

    TEST_FUNCTION_END(state);
}
