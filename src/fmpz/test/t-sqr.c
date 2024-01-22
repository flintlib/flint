/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_sqr, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        fmpz_t a, b;
        mpz_t c, d, e;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(c);
        mpz_init(d);
        mpz_init(e);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(c, a);

        aliasing = n_randint(state, 2);

        if (aliasing == 0)
        {
            fmpz_sqr(b, a);
        }
        else
        {
            fmpz_set(b, a);
            mpz_set(d, e);
            fmpz_sqr(b, b);
        }

        mpz_mul(d, c, c);

        fmpz_get_mpz(e, b);

        result = (mpz_cmp(d, e) == 0) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf(
                    "aliasing = %d\n"
                    "c = %Zd\n"
                    "expected: %Zd\n"
                    "got:      %Zd\n",
                    aliasing, c, d, e);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(c);
        mpz_clear(d);
        mpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
