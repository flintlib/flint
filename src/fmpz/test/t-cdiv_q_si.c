/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "fmpz.h"
#include "long_extras.h"

TEST_FUNCTION_START(fmpz_cdiv_q_si, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        slong b;
        fmpz_t a, c;
        mpz_t d, e, f, g;

        fmpz_init(a);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        b = z_randtest_not_zero(state);

        fmpz_get_mpz(d, a);
        flint_mpz_set_si(e, b);

        if (n_randint(state, 2))
        {
            fmpz_cdiv_q_si(c, a, b);
        }
        else /* test aliasing */
        {
            fmpz_set(c, a);
            fmpz_cdiv_q_si(c, c, b);
        }

        mpz_cdiv_q(f, d, e);

        fmpz_get_mpz(g, c);

        result = (mpz_cmp(f, g) == 0) && _fmpz_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL (1):\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd\n", d, e, f, g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}
