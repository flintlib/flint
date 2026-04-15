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

TEST_FUNCTION_START(n_mod_lemire, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong n, x, npre, r1, r2;

        do {
            n = n_randtest_not_zero(state);
        } while (n >= (UWORD(1) << (FLINT_BITS / 2)));

        npre = n_lemire_precomp(n);

        do {
            x = n_randtest_not_zero(state);
        } while (x >= (UWORD(1) << (FLINT_BITS / 2)));

        r1 = n_mod_lemire(x, n, npre);
        r2 = x % n;

        result = (r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = %wu, npre = %wu, x = %wu\n"
                    "r1 = %wu, r2 = %wu\n",
                    n, npre, x, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
