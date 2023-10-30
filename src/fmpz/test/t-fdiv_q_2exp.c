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

TEST_FUNCTION_START(fmpz_fdiv_q_2exp, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f;
        ulong x;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randint(state, 200);

        if (n_randint(state, 2))
        {
            fmpz_fdiv_q_2exp(b, a, x);
        }
        else /* test aliasing */
        {
            fmpz_set(b, a);
            fmpz_fdiv_q_2exp(b, b, x);
        }

        mpz_fdiv_q_2exp(e, d, x);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, exp = %Mu\n", d, e, f, x);
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
