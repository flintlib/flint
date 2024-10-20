/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "n_fft.h"

TEST_FUNCTION_START(n_fft_ctx_init2, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        // take random prime in [8, 2**61)
        ulong bits = 4 + n_randint(state, 57);
        ulong p = n_randprime(state, bits, 1);
        ulong max_depth = flint_ctz(p-1);

        // we need p such that 8 divides p-1
        while (max_depth < 3)
        {
            p = n_randprime(state, bits, 1);
            max_depth = flint_ctz(p-1);
        }

        // take depth between 0 and min(12, max_depth)
        ulong depth = n_randint(state, FLINT_MIN(10, max_depth));

        // init
        n_fft_ctx_t F;
        n_fft_ctx_init2(F, depth, p);

        if (F->mod != p)
            TEST_FUNCTION_FAIL(
                    "mod = %wu\n"
                    "F->mod = %wu\n",
                    p, F->mod);

        if (F->mod2 != 2*p)
            TEST_FUNCTION_FAIL(
                    "F->mod = %wu\n"
                    "F->mod2 = %wu\n",
                    F->mod, F->mod2);

        if (F->mod4 != 4*p)
            TEST_FUNCTION_FAIL(
                    "F->mod = %wu\n"
                    "F->mod4 = %wu\n",
                    F->mod, F->mod4);

        n_fft_ctx_clear(F);
    }

    TEST_FUNCTION_END(state);
}
