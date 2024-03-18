/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

TEST_FUNCTION_START(flint_clz, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n;
        unsigned int count = 0;

        n = n_randtest(state);

        if (n != 0)
            count = flint_clz(n);

        result = ((n == UWORD(0)) || (((slong)(n << count) < WORD(0)) && (r_shift(n, FLINT_BITS-count) == UWORD(0))));
        if (!result)
            TEST_FUNCTION_FAIL("n = %wu, count = %u\n", n, count);
    }

    TEST_FUNCTION_END(state);
}
