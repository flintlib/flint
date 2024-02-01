/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "gmpcompat.h"

TEST_FUNCTION_START(z_sizeinbase, state)
{
    int i, result;

    for (i = 0; i < 100000; i++)
    {
        slong a;
        mpz_t b;
        int base;
        size_t r1, r2;

        a = z_randtest(state);
        flint_mpz_init_set_si(b, a);
        base = (int) (n_randint(state, 61) + 2);

        r1 = z_sizeinbase(a, base);
        r2 = mpz_sizeinbase(b, base);
        result = (r1 == r2 || r1 + 1 == r2);

        if (!result)
            TEST_FUNCTION_FAIL(
                    "b = %{mpz}\n"
                    "base = %d\n"
                    "r1 = %wu\n"
                    "r2 = %wu\n",
                    b, base, (ulong) r1, (ulong) r2);

        mpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
