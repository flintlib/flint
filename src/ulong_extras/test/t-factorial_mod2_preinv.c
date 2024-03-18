/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

#ifndef n_factorial_mod2_foolproof
#define n_factorial_mod2_foolproof n_factorial_mod2_foolproof
static ulong
n_factorial_mod2_foolproof(ulong n, ulong p, ulong pinv)
{
    ulong prod = UWORD(1) % p;

    while (n)
    {
        prod = n_mulmod2_preinv(prod, n, p, pinv);
        n--;
    }

    return prod;
}
#endif

TEST_FUNCTION_START(n_factorial_mod2_preinv, state)
{
    slong ix;

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
        ulong n, p, pinv, x, y;
        int type;

        type = n_randint(state, 100);

        if (type == 0)
        {
            /* n is big */
            n = n_randint(state, 1000000);
        }
        else
        {
            /* n is small */
            n = n_randint(state, 1000);
        }

        p = n_randtest_not_zero(state);
        pinv = n_preinvert_limb(p);

        x = n_factorial_mod2_preinv(n, p, pinv);
        y = n_factorial_mod2_foolproof(n, p, pinv);

        if (x != y)
            TEST_FUNCTION_FAIL(
                    "n = %wu\n"
                    "p = %wu\n"
                    "x = %wu\n"
                    "y = %wu\n",
                    n, p, x, y);
    }

    TEST_FUNCTION_END(state);
}
