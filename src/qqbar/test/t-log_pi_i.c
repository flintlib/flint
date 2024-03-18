/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_log_pi_i, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        slong p, p2;
        ulong q, q2;
        int res;

        qqbar_init(x);
        qqbar_init(y);

        q = 1 + n_randint(state, 30);
        p = n_randint(state, 1000);
        p -= 500;

        qqbar_exp_pi_i(x, p, q);
        res = qqbar_log_pi_i(&p2, &q2, x);
        if (res)
            qqbar_exp_pi_i(y, p2, q2);

        if (res == 0 || !qqbar_equal(x, y) || n_gcd(FLINT_ABS(p2), q2) != 1 || p2 > (slong) q2 || p2 <= -(slong) q2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("res = %d\n\n", res);
            flint_printf("p, p2 = %wd %wd\n\n", p, p2);
            flint_printf("q, q2 = %wu %wu\n\n", q, q2);
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
    }

    TEST_FUNCTION_END(state);
}
