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
    int n;

    FLINT_TEST_INIT(state);
    
    flint_printf("prime_pi....");
    fflush(stdout);

    for (n=1; n<10000 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        if ((n_prime_pi(n-1)+1 == n_prime_pi(n)) != n_is_prime(n))
        {
            flint_printf("FAIL:\n");
            flint_printf("expected pi(%d) + 1 = pi(%d)\n", n-1, n); 
            fflush(stdout);
            flint_abort();
        }
    }

    for (n=1; n<5000 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        if (n_prime_pi(n_nth_prime(n)) != n)
        {
            flint_printf("FAIL:\n");
            flint_printf("expected pi(prime(%d)) = %d\n", n, n); 
            fflush(stdout);
            flint_abort();
        }
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
