/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "fmpz.h"
#include "aprcl.h"

int main(void)
{
    slong bits;
    slong allbits[] = {80, 96, 112, 128, 160, 192, 224, 256, 320, 384, 448, 512, 0};
    slong bi;

    flint_printf("per-proof time, same primes:\n");
    flint_printf("%6s %14s %14s %10s\n", "bits", "fmpz_is_prime", "aprcl_jacobi", "overhead");

    for (bi = 0; (bits = allbits[bi]) != 0; bi++)
    {
        fmpz_t n;
        double t1 = 1e300, t2 = 1e300;
        slong i, tries;

        fmpz_init(n);
        fmpz_one(n);
        fmpz_mul_2exp(n, n, bits - 1);
        fmpz_add_ui(n, n, 12345);
        fmpz_nextprime(n, n, 0);

        for (tries = 0; tries < 3; tries++)
        {
            timeit_t timer;
            slong reps = FLINT_MAX(1, 400 / FLINT_MAX(1, (bits * bits * bits) / 300000));

            timeit_start(timer);
            for (i = 0; i < reps; i++)
                if (fmpz_is_prime(n) != 1)
                    flint_abort();
            timeit_stop(timer);
            t1 = FLINT_MIN(t1, (double) timer->wall / reps);

            timeit_start(timer);
            for (i = 0; i < reps; i++)
                if (aprcl_is_prime_jacobi(n) != 1)
                    flint_abort();
            timeit_stop(timer);
            t2 = FLINT_MIN(t2, (double) timer->wall / reps);
        }

        flint_printf("%6wd %12.2f ms %12.2f ms %9.0f%%\n",
                     bits, t1, t2, 100.0 * (t1 - t2) / t2);

        fmpz_clear(n);
    }

    return 0;
}
