/*
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

TEST_FUNCTION_START(fmpz_remove, state)
{
    int i, result;

    /* Compare with GMP, random input */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, g;
        slong x, y;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest_not_zero(a, state, 200);
        do {
            fmpz_randtest_not_zero(b, state, 200);
            fmpz_abs(b, b);
        } while (fmpz_is_one(b));

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        aliasing = n_randint(state, 4);

        if (aliasing == 0)
        {
            x = fmpz_remove(c, a, b);
        }
        else if (aliasing == 1)
        {
            fmpz_set(a, b);
            mpz_set(d, e);
            x = fmpz_remove(c, a, a);
        }
        else if (aliasing == 2)
        {
            fmpz_set(c, a);
            x = fmpz_remove(c, c, b);
        }
        else
        {
            fmpz_set(c, b);
            x = fmpz_remove(c, a, c);
        }

        y = mpz_remove(f, d, e);

        fmpz_get_mpz(g, c);

        result = ((x == y) && (mpz_cmp(f, g) == 0)) && _fmpz_is_canonical(c);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Compare with GMP, random input but ensure that factors exist */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, pow;
        mpz_t d, e, f, g;
        slong x, y;
        ulong n;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(pow);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest_not_zero(a, state, 200);
        do {
            fmpz_randtest_not_zero(b, state, 200);
            fmpz_abs(b, b);
        } while (fmpz_is_one(b));

        n = n_randint(state, 10);
        fmpz_pow_ui(pow, b, n);
        fmpz_mul(a, a, pow);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        x = fmpz_remove(c, a, b);
        y = mpz_remove(f, d, e);

        fmpz_get_mpz(g, c);

        result = ((x == y) && (x >= n) && (mpz_cmp(f, g) == 0));

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(pow);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}
