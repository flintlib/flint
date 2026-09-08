/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(flint_mpn_divisible, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        fmpz_t a, b, t;
        mpz_t ma, mb;
        int res, res2;
        slong abits, bbits;

        fmpz_init(a); fmpz_init(b); fmpz_init(t); mpz_init(ma); mpz_init(mb);

        switch (n_randint(state, 4))
        {
            case 0: bbits = 1 + n_randint(state, 130); abits = 1 + n_randint(state, 260); break;
            case 1: bbits = 1 + n_randint(state, 2000); abits = 1 + n_randint(state, 4000); break;
            case 2: bbits = 1 + n_randint(state, 20000); abits = 1 + n_randint(state, 60000); break;
            default: bbits = 1 + n_randint(state, 300); abits = 1 + n_randint(state, 20000); break;
        }

        fmpz_randbits(b, state, bbits);
        fmpz_abs(b, b);
        if (fmpz_is_zero(b))
            fmpz_one(b);
        if (n_randint(state, 3) == 0)
            fmpz_mul_2exp(b, b, n_randint(state, 200));   /* even divisors, zero limbs */
        if (n_randint(state, 5) == 0)
            fmpz_mul_ui(b, b, 3 * 7 * 11);

        fmpz_randbits(a, state, abits);
        fmpz_abs(a, a);
        if (n_randint(state, 2))
        {
            /* divisible or nearly */
            fmpz_mul(a, a, b);
            if (n_randint(state, 3) == 0)
            {
                fmpz_randbits(t, state, n_randint(state, bbits + 1));
                fmpz_add(a, a, t);
            }
        }

        fmpz_get_mpz(ma, a);
        fmpz_get_mpz(mb, b);

        res2 = mpz_divisible_p(ma, mb);
        res = flint_mpn_divisible(ma->_mp_d, ma->_mp_size, mb->_mp_d, mb->_mp_size);

        if (res != res2)
            TEST_FUNCTION_FAIL("res = %d, res2 = %d\na = %{fmpz}\nb = %{fmpz}\n", res, res2, a, b);

        /* fmpz level */
        if (fmpz_divisible(a, b) != res2)
            TEST_FUNCTION_FAIL("fmpz_divisible: a = %{fmpz}\nb = %{fmpz}\n", a, b);

        fmpz_clear(a); fmpz_clear(b); fmpz_clear(t); mpz_clear(ma); mpz_clear(mb);
    }

    TEST_FUNCTION_END(state);
}
