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

TEST_FUNCTION_START(n_discrete_log_bsgs, state)
{
    int i;

    for (i = 0; i < 1000; i++)
    {
        mp_limb_t p, root, b, d, result;
        double pinv;

        p = n_randprime(state, 26, 1);
        pinv = n_precompute_inverse(p);
        root = n_primitive_root_prime(p);
        b = n_randint(state, p - 1) + 1;
        d = n_discrete_log_bsgs(b, root, p);

        result = n_powmod_precomp(root, d, p, pinv);

        if (result != b)
        {
            flint_printf("FAIL:\n");
            flint_printf("%wu ** (%wu) == %wu != %wu\n", root, d, result, b);
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
