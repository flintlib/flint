/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_mod_barrett, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong n, x, npre, r1, r2, r3;

        n = n_randtest_not_zero(state);
        if (n_randint(state, 8) == 0)
            n = UWORD(1) << n_randint(state, FLINT_BITS);

        if (n_randint(state, 2))
        {
            npre = (n == 1) ? UWORD_MAX : n_barrett_precomp(n);
        }
        else
        {
            /* floor((2^FLINT_BITS - 1) / n): one less than
               n_barrett_precomp(n) when n is a power of two */
            ulong npre2 = (n == 1) ? UWORD_MAX : n_barrett_precomp(n) - ((n & (n - 1)) == 0);

            npre = n_barrett_precomp_preinv(n, n_preinvert_limb(n));

            if (npre != npre2)
                TEST_FUNCTION_FAIL("n_barrett_precomp_preinv\n"
                    "n = %wu, npre = %wu, expected %wu\n", n, npre, npre2);
        }

        x = n_randtest(state);

        r1 = n_mod_barrett(x, n, npre);
        r2 = n_mod_barrett_lazy(x, n, npre);
        r3 = x % n;

        result = (r1 == r3) && (r2 == r3 || (r3 <= UWORD_MAX - n && r2 == r3 + n));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = %wu, npre = %wu, x = %wu\n"
                    "r1 = %wu, r2 = %wu, r3 = %wu\n",
                    n, npre, x, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
