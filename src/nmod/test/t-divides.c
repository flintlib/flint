/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"

TEST_FUNCTION_START(nmod_divides, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        mp_limb_t n, x, y, xy, z;
        int div;

        n = n_randtest_not_zero(state);

        nmod_init(&mod, n);

        x = n_randtest(state) % n;  
        y = n_randtest(state) % n;

        xy = nmod_mul(x, y, mod);

        div = nmod_divides(&z, xy, x, mod);

        /* Claimed divisible, so check this. */
        if (!div || nmod_mul(z, x, mod) != xy)
            TEST_FUNCTION_FAIL("n = %wu, div = %d, x = %wu, y = %wu, z = %wu\n", n, div, x, y, z);

        div = nmod_divides(&z, x, y, mod);

        /* If claimed not divisible, verify by brute force. */
        if (!div && n <= 50)
        {
            for (z = 0; z < n; z++)
            {
                if (nmod_mul(z, y, mod) == x)
                    TEST_FUNCTION_FAIL("n = %wu, div = %d, x = %wu, y = %wu, z = %wu\n", n, div, x, y, z);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
