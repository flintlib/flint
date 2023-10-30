/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_cdiv_ui, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t b;
        ulong x, r1, r2;

        fmpz_init(a);
        mpz_init(b);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(b, a);
        x = n_randtest_not_zero(state);

        r1 = fmpz_cdiv_ui(a, x);
        r2 = flint_mpz_cdiv_ui(b, x);

        result = (r1 == r2) && _fmpz_is_canonical(a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("b = %Zd, x = %wu, r1 = %wu, r2 = %wu\n", b, x, r1, r2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        mpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
