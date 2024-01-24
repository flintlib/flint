/*
    Copyright (C) 2009, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_lll_mod_preinv, state)
{
    int i, result;

    /* test n_lll_mod_preinv against n_ll_mod_preinv */
    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, nh, nm, nl, r1, r2, rm;

        d = n_randtest_not_zero(state);
        nh = n_randtest(state) % d;
        nm = n_randtest(state);
        nl = n_randtest(state);

        dinv = n_preinvert_limb(d);

        r2 = n_lll_mod_preinv(nh, nm, nl, d, dinv);

        rm = n_ll_mod_preinv(nh, nm, d, dinv);
        r1 = n_ll_mod_preinv(rm, nl, d, dinv);

        result = (r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "nh = %wu, nm = %wd, nl = %wu, d = %wu, dinv = %wu\n"
                    "r1 = %wu, r2 = %wu\n",
                    nh, nm, nl, d, dinv, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
