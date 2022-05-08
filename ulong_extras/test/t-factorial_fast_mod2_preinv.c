/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"


static ulong
n_factorial_mod2_foolproof(ulong n, ulong p, ulong pinv)
{
    ulong prod = UWORD(1) % p;

    while (n)
    {
        prod = n_mulmod2_preinv(prod, n, p, pinv);
        n--;
    }

    return prod;
}

int main(void)
{
    ulong n;
    int j;

    FLINT_TEST_INIT(state);    

    flint_printf("factorial_fast_mod2_preinv....");
    fflush(stdout);

    for (n = 0; n < 100 * flint_test_multiplier(); n++)
    {
        ulong p, pinv, x, y;

        for (j = 0; j < 10; j++)
        {
            p = n_randtest_not_zero(state);
            pinv = n_preinvert_limb(p);
            x = n_factorial_fast_mod2_preinv(n, p, pinv);
            y = n_factorial_mod2_foolproof(n, p, pinv);

            if (x != y)
            {
                flint_printf("FAIL:\n");
                flint_printf("n = %wu\np = %wu\nx = %wu\ny = %wu\n", n, p, x, y);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
