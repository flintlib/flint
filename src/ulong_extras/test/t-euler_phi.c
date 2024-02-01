/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_euler_phi, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        ulong nx, kx;
        ulong r1, r2;

        nx = n_randbits(state, n_randint(state, 18));
        r1 = 0;

        for (kx = 1; kx <= nx; kx++)
            r1 += (n_gcd(nx, kx) == 1);

        r2 = n_euler_phi(nx);

        result = (r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "phi(%wu) = %wu, got %wu\n",
                    nx, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
