/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_primitive_root_prime, state)
{
    int i, j;

    for (i = 0; i < 100; i++)
    {
        n_factor_t factors;
        mp_limb_t p, root;
        double pinv;

        n_factor_init(&factors);
        p = n_randtest_prime(state, 1);
        pinv = n_precompute_inverse(p);
        n_factor(&factors, p - 1, 1);

        root = n_primitive_root_prime(p);

        for (j = 0; j < factors.num; j++)
        {
            if (n_powmod_precomp(root, (p-1) / factors.p[j], p, pinv) == 1)
            {
                flint_printf("FAIL:\n");
                flint_printf("%wu ** (%wu / %wu) == 1 mod %wu\n", root, p-1, factors.p[j], p);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    TEST_FUNCTION_END(state);
}
