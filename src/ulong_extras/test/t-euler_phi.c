/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int n, k, t;

    FLINT_TEST_INIT(state);

    flint_printf("euler_phi....");
    fflush(stdout);

    for (n = 0; n < 20 * flint_test_multiplier(); n++)
    {
        t = 0;
        for (k = 1; k <= n; k++)
            t += (n_gcd(n, k) == 1);
        if (t != n_euler_phi(n))
        {
            flint_printf("FAIL:\n");
            flint_printf("phi(%d) = %d, got %wu\n", n, t, n_euler_phi(n)); 
            fflush(stdout);
            flint_abort();
        }
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
