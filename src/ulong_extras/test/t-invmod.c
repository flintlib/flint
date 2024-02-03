/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_invmod, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        ulong a, b, t, r, binv, ph, pl;

        do
        {
            a = n_randtest(state);
            b = n_randtest(state);
        } while ((a >= b) || (n_gcd(b, a) != UWORD(1)));

        t = n_invmod(a, b);

        binv = n_preinvert_limb(b);
        umul_ppmm(ph, pl, t, a);
        r = n_ll_mod_preinv(ph, pl, b, binv);

        result = ((r == UWORD(0) && b == UWORD(1)) || r == UWORD(1));
        if (!result)
            TEST_FUNCTION_FAIL("a = %wu, b = %wu, r = %wd\n", a, b, r);
    }

    TEST_FUNCTION_END(state);
}
