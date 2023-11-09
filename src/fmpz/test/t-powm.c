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
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_powm, state)
{
    int i, result;

    /* Compare with GMP */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, m;
        fmpz_t x;
        mpz_t y;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(x);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(m);
        mpz_init(y);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(c, state, 200);
        fmpz_abs(c, c);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(m, c);
        fmpz_randtest_unsigned(x, state, 20);
        fmpz_get_mpz(y, x);

        aliasing = n_randint(state, 5);

        if (aliasing == 0)
        {
            fmpz_powm(b, a, x, c);
        }
        else if (aliasing == 1)
        {
            fmpz_set(b, a);
            fmpz_powm(b, a, x, c);
        }
        else if (aliasing == 2)
        {
            fmpz_set(b, c);
            fmpz_powm(b, a, x, b);
        }
        else if (aliasing == 3)
        {
            fmpz_set(a, c);
            mpz_set(d, m);
            fmpz_powm(b, a, x, a);
        }
        else
        {
            fmpz_set(b, c);
            mpz_set(d, m);
            fmpz_powm(b, b, x, b);
        }

        mpz_powm(e, d, y, m);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL (cmp f with GMP e := d^y mod m):\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, y = %Zd, m = %Zd\n", d, e, f, y, m);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(x);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(m);
        mpz_clear(y);
    }

    TEST_FUNCTION_END(state);
}
