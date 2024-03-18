/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "fmpz.h"
#include "long_extras.h"

TEST_FUNCTION_START(fmpz_mul_si, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f;
        slong x;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = z_randtest(state);

        if (n_randint(state, 2))
        {
            fmpz_mul_si(b, a, x);
        }
        else
        {
            fmpz_set(b, a);
            fmpz_mul_si(b, b, x);
        }

        flint_mpz_mul_si(e, d, x);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, x = %Md\n", d, e, f, x);
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
