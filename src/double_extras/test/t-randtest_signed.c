/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <float.h>
#include "ulong_extras.h"
#include "double_extras.h"

TEST_FUNCTION_START(d_randtest_signed, state)
{
    double x;
    slong iter;

    /* check that values lie in [0.5, 1) U {0} for minexp = maxexp = 0 */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        x = d_randtest_signed(state, 0, 0);
        if ((fabs(x) < 0.5 && x != 0) || fabs(x) >= 1)
            TEST_FUNCTION_FAIL("x = %.17g\n", x);
    }

    TEST_FUNCTION_END(state);
}
