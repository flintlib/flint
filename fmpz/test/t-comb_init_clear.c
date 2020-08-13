/*
    Copyright (C) 2008, 2009, William Hart
    Copyright (C) 2010 Fredrik Johansson

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
#include "fmpz.h"


int main()
{
    slong i, j;
    mp_limb_t n;
    slong num_primes;
    mp_limb_t * primes;
    mp_limb_t p;
    fmpz_comb_t comb;
    FLINT_TEST_INIT(state);
    

    flint_printf("comb_init/clear....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        n = n_randint(state, 10);
        num_primes = (WORD(1) << n);
        primes = (mp_limb_t *) flint_malloc(num_primes * sizeof(mp_limb_t));
        p = n_nextprime((UWORD(1) << (FLINT_BITS-1)) - WORD(10000000), 0);

        for (j = 0; j < num_primes; j++)
        {
            primes[j] = p;
            p = n_nextprime(p, 0);
        }

        fmpz_comb_init(comb, primes, num_primes);
        fmpz_comb_clear(comb);
        flint_free(primes);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
