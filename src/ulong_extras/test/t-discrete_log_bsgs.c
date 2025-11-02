/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

#define big_prime(state) n_randprime(state, 22 + n_randint(state, 8), 0)
#define small_prime(state) n_randprime(state, 2 + n_randint(state, 18), 0)

TEST_FUNCTION_START(n_discrete_log_bsgs, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p, root, b, b2, d;
        double pinv;

        p = n_randint(state, 128) ? small_prime(state) : big_prime(state);
        pinv = n_precompute_inverse(p);
        root = n_primitive_root_prime(p);
        b = n_randint(state, p - 1) + 1;
        d = n_discrete_log_bsgs(b, root, p);
        b2 = n_powmod_precomp(root, d, p, pinv);

        result = b == b2;

        if (!result)
            TEST_FUNCTION_FAIL("%wu ** (%wu) == %wu != %wu\n", root, d, b2, b);
    }

    TEST_FUNCTION_END(state);
}

#undef big_prime
#undef small_prime
