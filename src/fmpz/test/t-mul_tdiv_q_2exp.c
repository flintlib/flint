/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_mul_tdiv_q_2exp, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, g;
        ulong exp;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);
        exp = n_randint(state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        aliasing = n_randint(state, 4);

        if (aliasing == 0)
        {
            fmpz_mul_tdiv_q_2exp(c, a, b, exp);
        }
        else if (aliasing == 1)
        {
            fmpz_set(a, b);
            mpz_set(d, e);
            fmpz_mul_tdiv_q_2exp(c, a, a, exp);
        }
        else if (aliasing == 2)
        {
            fmpz_set(c, a);
            fmpz_mul_tdiv_q_2exp(c, c, b, exp);
        }
        else
        {
            fmpz_set(c, b);
            fmpz_mul_tdiv_q_2exp(c, a, c, exp);
        }

        mpz_mul(f, d, e);
        mpz_tdiv_q_2exp(f, f, exp);

        fmpz_get_mpz(g, c);

        result = (mpz_cmp(f, g) == 0) && _fmpz_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd, exp = %Mu\n",
                d, e, f, g, exp);
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

    TEST_FUNCTION_END(state);
}
