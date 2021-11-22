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

    flint_printf("sqrtmod....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random integers */
    {
        mp_limb_t a, b, p, pinv;

        p = n_randtest_prime(state, 0);
        a = n_randtest(state) % p;

        b = n_sqrtmod(a, p);
        pinv = n_preinvert_limb(p);

        result = (b == 0 || n_mulmod2_preinv(b, b, p, pinv) == a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wu\n", p);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            fflush(stdout);
            flint_abort();
        }
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares */
    {
        mp_limb_t a, b, p, pinv;

        p = n_randtest_prime(state, 0);

        do 
            b = n_randtest(state) % p;
        while (b == 0);

        pinv = n_preinvert_limb(p);
        a = n_mulmod2_preinv(b, b, p, pinv);

        b = n_sqrtmod(a, p);

        result = (n_mulmod2_preinv(b, b, p, pinv) == a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wu\n", p);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            fflush(stdout);
            flint_abort();
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
