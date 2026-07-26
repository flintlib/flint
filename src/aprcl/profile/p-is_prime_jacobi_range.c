/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "fmpz.h"
#include "aprcl.h"

int main(int argc, char * argv[])
{
    slong bits;
    flint_rand_t state;

    flint_rand_init(state);
    flint_rand_set_seed(state, 12345, 67890);

    flint_printf("aprcl_is_prime_jacobi, best of 3, per-proof time:\n");

    for (bits = 80; bits <= 512; bits += (bits < 128) ? 16 : ((bits < 256) ? 32 : 64))
    {
        fmpz_t n;
        slong i;
        double t, best = 1e100;

        fmpz_init(n);

        /* deterministic prime of the given size */
        fmpz_one(n);
        fmpz_mul_2exp(n, n, bits - 1);
        fmpz_add_ui(n, n, 12345);
        fmpz_nextprime(n, n, 0);

        for (i = 0; i < 3; i++)
        {
            timeit_t timer;
            slong reps, r;
            reps = FLINT_MAX(1, 400 / FLINT_MAX(1, (bits * bits * bits) / 300000));
            timeit_start(timer);
            for (r = 0; r < reps; r++)
                if (!aprcl_is_prime_jacobi(n))
                    flint_abort();
            timeit_stop(timer);
            t = (double) timer->wall / reps;
            if (t < best)
                best = t;
        }

        flint_printf("bits = %3wd:  %8.2f ms\n", bits, best);

        fmpz_clear(n);
    }

    flint_rand_clear(state);
    return 0;
}
