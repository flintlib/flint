/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(compute_primes, state)
{
    slong i, lim = 1000000;
    n_primes_t pg;
    ulong * ref_primes;
    double * ref_inverses;

    ref_primes = flint_malloc(sizeof(ulong) * lim);
    ref_inverses = flint_malloc(sizeof(double) * lim);

    n_primes_init(pg);
    for (i = 0; i < lim; i++)
    {
        ref_primes[i] = n_primes_next(pg);
        ref_inverses[i] = n_precompute_inverse(ref_primes[i]);
    }
    n_primes_clear(pg);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong n;
        const ulong * primes;
        const double * inverses;

        n = n_randint(state, lim);

        primes = n_primes_arr_readonly(n + 1);
        inverses = n_prime_inverses_arr_readonly(n + 1);

        if (primes[n] != ref_primes[n] || inverses[n] != ref_inverses[n])
            TEST_FUNCTION_FAIL(
                    "n = %wd, p1 = %wu, p2 = %wu\n"
                    "inv1 = %g, inv2 = %g\n",
                    n, primes[n], ref_primes[n], inverses[n], ref_inverses[n]);

        if (n_randint(state, 50) == 0)
            n_cleanup_primes();
    }

    flint_free(ref_primes);
    flint_free(ref_inverses);

    TEST_FUNCTION_END(state);
}
