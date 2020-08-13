/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("clog....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t a = 0, b = 0, k, x;

        while (a < 1)
            a = n_randtest(state);
        while (b < 2)
            b = n_randtest(state);

        k = n_clog(a, b);
        x = (k > 0) ? n_pow(b, k - 1) : 0;

        result = (x < a && ((k > 0) ? ((a + b - 1) / b <= x) : (a <= 1)));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            flint_printf("x = %wu\n", x);
            flint_printf("k = %wu\n", k);
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
