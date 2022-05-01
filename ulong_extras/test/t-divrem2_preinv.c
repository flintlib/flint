/*
    Copyright (C) 2009, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("divrem2_preinv....");
    fflush(stdout);

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, n, q1, q2, r1, r2;

        d = n_randtest_not_zero(state);
        n = n_randtest(state);

        dinv = n_preinvert_limb(d);

        r1 = n_divrem2_preinv(&q1, n, d, dinv);
        q2 = n / d;
        r2 = n % d;

        result = (q1 == q2 && r1 == r2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, d = %wu, dinv = %wu\n", n, d, dinv);
            flint_printf("q1 = %wu, q2 = %wu\n", q1, q2);
            flint_printf("r1 = %wu, r2 = %wu\n", r1, r2);
            fflush(stdout);
            flint_abort();
        }
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
