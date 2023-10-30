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

/* Use the definition of GMP versions >= 6.0 */
int
mpz_invert2(mpz_t a, const mpz_t b, const mpz_t c)
{
    if (mpz_cmpabs_ui(c, 1) == 0)
    {
        mpz_set_ui(a, 0);
        return 1;
    }
    else
        return mpz_invert(a, b, c);
}

TEST_FUNCTION_START(fmpz_invmod, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, g;
        int r1, r2;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        aliasing = n_randint(state, 4);

        if (aliasing == 0)
        {
            r1 = fmpz_invmod(c, a, b);
        }
        else if (aliasing == 1)
        {
            fmpz_set(a, b);
            mpz_set(d, e);
            r1 = fmpz_invmod(c, a, a);
        }
        else if (aliasing == 2)
        {
            fmpz_set(c, a);
            r1 = fmpz_invmod(c, c, b);
        }
        else
        {
            fmpz_set(c, b);
            r1 = fmpz_invmod(c, a, c);
        }

        r2 = mpz_invert2(f, d, e);

        fmpz_get_mpz(g, c);

        result = (r1 != 0 && r2 != 0 && (mpz_cmp(f, g) == 0)) || (r1 == 0 && r2 == 0);
        result = result && _fmpz_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, r1 = %d, r2 = %d\n", d,
                 e, f, g, r1, r2);
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
