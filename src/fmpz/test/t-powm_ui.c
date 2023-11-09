/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_powm_ui, state)
{
    int i, result;

    /* Compare with GMP */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, m;
        ulong x;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(m);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(m, c);
        x = n_randtest(state);

        aliasing = n_randint(state, 4);

        if (aliasing == 0)
        {
            fmpz_powm_ui(b, a, x, c);
        }
        else if (aliasing == 1)
        {
            fmpz_set(b, a);
            fmpz_powm_ui(b, b, x, c);
        }
        else if (aliasing == 2)
        {
            fmpz_set(a, c);
            mpz_set(d, m);
            fmpz_powm_ui(b, a, x, a);
        }
        else if (aliasing == 3)
        {
            fmpz_set(b, c);
            mpz_set(d, m);
            fmpz_powm_ui(b, b, x, b);
        }

        flint_mpz_powm_ui(e, d, x, m);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, x = %Mu, m = %Zd\n", d, e, f, x, m);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(m);
    }

    TEST_FUNCTION_END(state);
}
