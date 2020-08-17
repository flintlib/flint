/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "ulong_extras.h"
#include "double_extras.h"

int
main(void)
{
    double x;
    slong iter;

    FLINT_TEST_INIT(state);

    flint_printf("randtest_signed....");
    fflush(stdout);

    /* check that values lie in [0.5, 1) U {0} for minexp = maxexp = 0 */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        x = d_randtest_signed(state, 0, 0);
        if ((fabs(x) < 0.5 && x != 0) || fabs(x) >= 1)
        {
            flint_printf("FAIL\n");
            flint_printf("x = %.17g\n", x);
            abort();
        }
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
