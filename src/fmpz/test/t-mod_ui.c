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

TEST_FUNCTION_START(fmpz_mod_ui, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f;
        ulong x, r1, r2;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randtest_not_zero(state);

        if (n_randint(state, 2))
        {
            r1 = fmpz_mod_ui(b, a, x);
        }
        else
        {
            fmpz_set(b, a);
            r1 = fmpz_mod_ui(b, b, x);
        }

        r2 = flint_mpz_fdiv_r_ui(e, d, x);

        fmpz_get_mpz(f, b);

        result = ((mpz_cmp(e, f) == 0) && (r1 == r2)) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, x = %wu, r1 = %wu, r2 = %wu\n", d,
                 e, f, x, r1, r2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    TEST_FUNCTION_END(state);
}
