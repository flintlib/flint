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

TEST_FUNCTION_START(fmpz_submul)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, g;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);
        fmpz_randtest(c, state, 200);

        if (n_randint(state, 10) == 0)
            fmpz_addmul(c, a, b);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);
        fmpz_get_mpz(f, c);

        fmpz_submul(c, a, b);
        mpz_submul(f, d, e);

        fmpz_get_mpz(g, c);

        result = (mpz_cmp(f, g) == 0);
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

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        mpz_t d, f, g;

        fmpz_init(a);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(c, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(f, c);

        fmpz_submul(c, a, a);
        mpz_submul(f, d, d);

        fmpz_get_mpz(g, c);

        result = (mpz_cmp(f, g) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, f = %Zd, g = %Zd\n", d, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Test aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_submul(a, a, b);
        mpz_submul(d, d, e);

        fmpz_get_mpz(f, a);

        result = (mpz_cmp(d, f) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd\n", d, e, f);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    /* Test aliasing of b and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_submul(b, a, b);
        mpz_submul(e, d, e);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(f, e) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd\n", d, e, f);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    TEST_FUNCTION_END;
}
