/*
    Copyright (C) 2012 Fredrik Johansson

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
    slong j, k, l;

    FLINT_TEST_INIT(state);   

    flint_printf("primes_jump_after....");
    fflush(stdout);

    for (j = 0; j < 10; j++)
    {
        n_primes_t iter;

        n_primes_init(iter);

        for (k = 0; k < 100; k++)
        {
            mp_limb_t p, q;

            q = n_randtest(state) % UWORD(1000000000);

            n_primes_jump_after(iter, q);

            for (l = 0; l < 100; l++)
            {
                p = n_primes_next(iter);
                q = n_nextprime(q, 0);

                if (p != q)
                {
                    flint_printf("FAIL\n");
                    flint_printf("p = %wu, q = %wu\n", p, q);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        n_primes_clear(iter);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
