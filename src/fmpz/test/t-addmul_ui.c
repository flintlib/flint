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

TEST_FUNCTION_START(fmpz_addmul_ui, state)
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
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);
        x = n_randtest(state);

        if (n_randint(state, 2))
        {
            fmpz_addmul_ui(b, a, x);
        }
        else  /* test aliasing */
        {
            fmpz_set(b, a);
            mpz_set(e, d);
            fmpz_addmul_ui(b, b, x);
        }

        flint_mpz_addmul_ui(e, d, x);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0) && _fmpz_is_canonical(b);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, x = %Mu\n", d, e, f, x);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t d, e;
        ulong x;

        fmpz_init(a);

        mpz_init(d);
        mpz_init(e);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randtest(state);

        fmpz_addmul_ui(a, a, x);
        flint_mpz_addmul_ui(d, d, x);

        fmpz_get_mpz(e, a);

        result = (mpz_cmp(d, e) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, x = %Mu\n", d, e, x);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);

        mpz_clear(d);
        mpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
