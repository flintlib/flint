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

TEST_FUNCTION_START(udiv_qrnnd, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, nh, nl, q, r, ph, pl;

        do
        {
            d = n_randtest_not_zero(state);
            nh = n_randtest(state);
        } while (nh >= d);
        nl = n_randtest(state);

        udiv_qrnnd(q, r, nh, nl, d);
        umul_ppmm(ph, pl, d, q);
        add_ssaaaa(ph, pl, ph, pl, UWORD(0), r);

        result = ((ph == nh) && (pl == nl));

        if (!result)
            TEST_FUNCTION_FAIL(
                    "nh = %wu, nl = %wu, d = %wu\n"
                    "ph = %wu, pl = %wu\n",
                    nh, nl, d, ph, pl);
    }

    TEST_FUNCTION_END(state);
}
