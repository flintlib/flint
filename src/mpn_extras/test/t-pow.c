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

TEST_FUNCTION_START(flint_mpn_pow, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_t x, y;
        mpz_t mx, my;
        mp_ptr r;
        mp_size_t xn, bound, rn, i;
        ulong e;
        slong xbits;

        fmpz_init(x); fmpz_init(y); mpz_init(mx); mpz_init(my);

        switch (n_randint(state, 4))
        {
            case 0: xbits = 1 + n_randint(state, 64); e = 1 + n_randint(state, 200); break;
            case 1: xbits = 1 + n_randint(state, 200); e = 1 + n_randint(state, 60); break;
            case 2: xbits = 1 + n_randint(state, 3000); e = 1 + n_randint(state, 12); break;
            default: xbits = 1 + n_randint(state, 8); e = n_randint(state, 5000); break;
        }

        fmpz_randbits(x, state, xbits);
        fmpz_abs(x, x);
        if (fmpz_is_zero(x))
            fmpz_one(x);
        if (n_randint(state, 4) == 0)
            fmpz_mul_2exp(x, x, n_randint(state, 150));    /* zero bits and limbs */
        if (n_randint(state, 20) == 0)
            fmpz_one(x);
        if (n_randint(state, 20) == 0)
            fmpz_set_ui(x, 2);

        fmpz_get_mpz(mx, x);
        xn = mx->_mp_size;

        bound = flint_mpn_pow_bound_limbs(mx->_mp_d, xn, e);
        r = flint_malloc((bound + 1) * sizeof(mp_limb_t));
        r[bound] = UWORD(0xdead);

        rn = flint_mpn_pow(r, mx->_mp_d, xn, e);
        mpz_pow_ui(my, mx, e);

        if (r[bound] != UWORD(0xdead))
            TEST_FUNCTION_FAIL("overrun: xbits = %wd, e = %wu, bound = %wd\n", xbits, e, bound);
        if (rn != my->_mp_size || rn > bound)
            TEST_FUNCTION_FAIL("size: rn = %wd, expected %wd, bound %wd\nx = %{fmpz}\ne = %wu\n", rn, (slong) my->_mp_size, bound, x, e);
        for (i = 0; i < rn; i++)
            if (r[i] != my->_mp_d[i])
                TEST_FUNCTION_FAIL("value: limb %wd\nx = %{fmpz}\ne = %wu\n", i, x, e);
        if (bound > rn + 3 && bound > 2 * rn && rn > 64)   /* the cheap bound for small results may be 2x loose */
            TEST_FUNCTION_FAIL("loose bound: rn = %wd, bound = %wd\nx = %{fmpz}\ne = %wu\n", rn, bound, x, e);

        flint_free(r);
        fmpz_clear(x); fmpz_clear(y); mpz_clear(mx); mpz_clear(my);
    }

    TEST_FUNCTION_END(state);
}
