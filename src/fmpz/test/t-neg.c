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

TEST_FUNCTION_START(fmpz_neg, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(c);
        mpz_init(d);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(c, a);

        if (n_randint(state, 2))
        {
            fmpz_neg(b, a);
        }
        else
        {
            fmpz_set(b, a);
            fmpz_neg(b, b);
        }

        mpz_neg(c, c);

        fmpz_get_mpz(d, b);

        result = (mpz_cmp(c, d) == 0) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(c);
        mpz_clear(d);
    }

    TEST_FUNCTION_END(state);
}
