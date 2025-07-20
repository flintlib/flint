/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_char_set_slong_vec, state)
{
    slong iter;

    /* Test: inverse of char_bit */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        ulong a, b, test, ta, tb;
        slong * vec;
        slong k;

        vec = flint_malloc(2 * g * sizeof(slong));
        a = n_randint(state, 1 << g);
        b = n_randint(state, 1 << g);

        for (k = 0; k < g; k++)
        {
            vec[k] = acb_theta_char_bit(a, k, g);
            vec[k + g] = acb_theta_char_bit(b, k, g);
        }

        test = acb_theta_char_set_slong_vec(vec, 2 * g);
        ta = acb_theta_char_set_slong_vec(vec, g);
        tb = acb_theta_char_set_slong_vec(vec + g, g);

        if (test != (a << g) + b || ta != a || tb != b)
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, test = %wd, a = %wd, b = %wd, ta = %wd, tb = %wd\n",
                g, test, a, b, ta, tb);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
