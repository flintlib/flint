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

TEST_FUNCTION_START(d_is_nan, state)
{
    double x;
    slong iter;

    /* check non-zero value returned if x == NaN */
    x = D_NAN;
    if (!d_is_nan(x))
    {
        flint_printf("FAIL\n");
        flint_printf("0 returned for %g\n", x);
        fflush(stdout);
        flint_abort();
    }

    /* check 0 returned if x != NaN */
    x = D_INF;
    if (d_is_nan(x))
    {
        flint_printf("FAIL\n");
        flint_printf("Non-zero returned for %g\n", x);
        fflush(stdout);
        flint_abort();
    }
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        x = d_randtest_signed(state, 0, 0);
        if (d_is_nan(x))
        {
            flint_printf("FAIL\n");
            flint_printf("Non-zero returned for %g\n", x);
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
