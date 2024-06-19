/*
    Copyright 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "ulong_extras.h"

slong bits[] = { 1, 32, 64, 128, 192, 224 };
slong reps[] = { 1000000, 100000, 10000, 100, 30, 10 };

int main(void)
{
    slong i, b;
    fmpz_t n;
    fmpz_factor_t fac;

    flint_set_num_threads(8);

    fmpz_init(n);
    fmpz_factor_init(fac);

    for (b = 0; b < 6; b++)
    {
        flint_printf("%wd x %wd bits:   ", reps[b], bits[b]);
        fmpz_one(n);
        fmpz_mul_2exp(n, n, bits[b] - 1);
        TIMEIT_ONCE_START
        for (i = 0; i < reps[b]; i++)
        {
            fmpz_add_ui(n, n, 1);
            fmpz_factor(fac, n);

/*
            if (bits[b] >= 190)
            {
                fmpz_factor_print(fac);
                flint_printf("\n");
            }
*/
        }
        TIMEIT_ONCE_STOP
    }

    fmpz_factor_clear(fac);
    fmpz_clear(n);

    flint_cleanup_master();
    return 0;
}
