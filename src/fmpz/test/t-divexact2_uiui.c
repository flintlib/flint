/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

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

TEST_FUNCTION_START(fmpz_divexact2_uiui, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        mpz_t e, f, g;
        ulong n, m;

        fmpz_init(a);
        fmpz_init(c);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        n = n_randtest_not_zero(state);
        m = n_randtest_not_zero(state);

        fmpz_mul_ui(c, a, n);
        fmpz_mul_ui(c, c, m);

        fmpz_get_mpz(e, c);

        if (n_randint(state, 2))
        {
            fmpz_divexact2_uiui(a, c, n, m);
        }
        else
        {
            fmpz_set(a, c);
            fmpz_divexact2_uiui(a, a, n, m);
        }

        flint_mpz_divexact_ui(f, e, n);
        flint_mpz_divexact_ui(f, f, m);

        fmpz_get_mpz(g, a);

        result = (mpz_cmp(f, g) == 0) && _fmpz_is_canonical(a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("n = %Mu, m = %Mu, e = %Zd, f = %Zd, g = %Zd\n",
                n, m, e, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}
