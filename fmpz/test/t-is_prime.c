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
    int i, result, r1;
    FLINT_TEST_INIT(state);

    flint_printf("is_prime....");
    fflush(stdout);

    /* test primes always pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, F;

        fmpz_init(p);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 160) + 2);
           fmpz_abs(p, p);
        } while (!fmpz_is_probabprime(p));

        r1 = fmpz_is_prime(p);

        result = (r1 == 1 || r1 == -1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_clear(F);
    }

    /* test composites never pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, a, F;

        fmpz_init(p);
        fmpz_init(a);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 80) + 2);
        } while (fmpz_cmp_ui(p, 2) < 0);
        do {
           fmpz_randbits(a, state, n_randint(state, 80) + 2);
        } while (fmpz_cmp_ui(a, 2) < 0);
        
        fmpz_mul(p, p, a);

        r1 = fmpz_is_prime(p);

        result = (r1 == 0 || r1 == -1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(F);
    }

    /* test issue 345 */
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set_ui(t, 13567);
        fmpz_pow_ui(t, t, 145);

        if (fmpz_is_prime(t) != 0)
        {
            flint_printf("FAIL:\n");
            fmpz_print(t); printf("\n");
            abort();
        }

        fmpz_clear(t);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
