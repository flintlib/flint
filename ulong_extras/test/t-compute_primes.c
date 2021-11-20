/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Fredrik Johansson

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

int main()
{
    slong i, lim = 1000000;
    n_primes_t pg;
    mp_limb_t * ref_primes;
    double * ref_inverses;
    FLINT_TEST_INIT(state);

    flint_printf("compute_primes....");
    fflush(stdout);
    

    ref_primes = flint_malloc(sizeof(mp_limb_t) * lim);
    ref_inverses = flint_malloc(sizeof(double) * lim);

    n_primes_init(pg);
    for (i = 0; i < lim; i++)
    {
        ref_primes[i] = n_primes_next(pg);
        ref_inverses[i] = n_precompute_inverse(ref_primes[i]);
    }
    n_primes_clear(pg);

    for (i = 0; i < 250; i++)
    {
        slong n;
        const mp_limb_t * primes;
        const double * inverses;

        n = n_randtest(state) % lim;

        primes = n_primes_arr_readonly(n + 1);
        inverses = n_prime_inverses_arr_readonly(n + 1);

        if (primes[n] != ref_primes[n] || inverses[n] != ref_inverses[n])
        {
            flint_printf("FAIL!\n");
            flint_printf("n = %wd, p1 = %wu, p2 = %wu\n", n, primes[n], ref_primes[n]);
            flint_printf("inv1 = %g, inv2 = %g\n", inverses[n], ref_inverses[n]);
            fflush(stdout);
            flint_abort();
        }

        if (n_randint(state, 20) == 0)
        {
            n_cleanup_primes();
        }
    }

    flint_free(ref_primes);
    flint_free(ref_inverses);
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

