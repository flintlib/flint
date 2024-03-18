/*
    Copyright 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "ulong_extras.h"

slong bits[] = { 1, 64, 128, 256, 1024, 4096 };
slong reps[] = { 10000000, 10000000, 1000000, 100000, 10000, 1000 };

int main(void)
{
    slong i, b, count;
    fmpz_t n;

    flint_set_num_threads(8);

    fmpz_init(n);

    for (b = 0; b < 6; b++)
    {
        count = 0;
        flint_printf("%wd x %wd bits:   ", reps[b], bits[b]);
        fmpz_one(n);
        fmpz_mul_2exp(n, n, bits[b] - 1);
        TIMEIT_ONCE_START
        for (i = 0; i < reps[b]; i++)
        {
            fmpz_add_ui(n, n, 1);
            count += fmpz_is_probabprime(n);
        }
        flint_printf("%wd found   ", count);
        TIMEIT_ONCE_STOP
    }

    fmpz_clear(n);

    flint_cleanup_master();
    return 0;
}
