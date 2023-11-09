/*
    Copyright (C) 2008, 2009, William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_comb_init_clear, state)
{
    slong i, j;
    mp_limb_t n;
    slong num_primes;
    mp_limb_t * primes;
    mp_limb_t p;
    fmpz_comb_t comb;

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

    TEST_FUNCTION_END(state);
}
