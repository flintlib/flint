/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result, count = 0;
    FLINT_TEST_INIT(state);

    flint_printf("is_strong_probabprime....");
    fflush(stdout);

    /* test primes always pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, b, F;

        fmpz_init(p);
        fmpz_init(b);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 330) + 2);
           fmpz_abs(p, p);
        } while (!fmpz_is_probabprime(p) || fmpz_cmp_ui(p, 2) == 0);

        do {
           fmpz_randbits(b, state, n_randint(state, 100) + 1);
           fmpz_abs(b, b);
        } while (fmpz_is_zero(b) || fmpz_is_one(b)); 
        
        result = fmpz_is_strong_probabprime(p, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            fmpz_print(b); printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_clear(b);
        fmpz_clear(F);
    }

    /* test composites rarely pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, a, b, F;

        fmpz_init(p);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 100) + 2);
        } while (fmpz_cmp_ui(p, 2) < 0);
        do {
           fmpz_randbits(a, state, n_randint(state, 100) + 2);
        } while (fmpz_cmp_ui(a, 2) < 0);
        
        do {
           fmpz_randbits(b, state, n_randint(state, 100) + 1);
           fmpz_abs(b, b);
        } while (fmpz_cmp_ui(b, 2) < 0); 
        
        fmpz_mul(p, p, a);

        if (fmpz_is_strong_probabprime(p, b))
           count++;

        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(F);
    }

    result = (count < flint_test_multiplier());
    if (!result)
    {
        flint_printf("FAIL:\n");
        flint_printf("count = %ld\n", count);
        abort();
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
