/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <float.h>
#include "ulong_extras.h"
#include "double_extras.h"

TEST_FUNCTION_START(d_log2, state)
{
    double x, res1, res2;
    slong iter;

    /* check change of base identity */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        x = d_randtest(state);
        res1 = d_log2(x) * log(2);
        res2 = log(x);
        if (fabs(res1 - res2) > D_EPS)
        {
            flint_printf("FAIL\n");
            flint_printf("x = %.20g\n", x);
            flint_printf("res1 = %.20g\n", res1);
            flint_printf("res2 = %.20g\n", res2);
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
