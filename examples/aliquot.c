/*
    Copyright (C) 2024 Matthias Gessinger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <flint/fmpz.h>
#include <flint/profiler.h>

int
main(int argc, char * argv[])
{
    fmpz_t n;
    ulong num_iterations = 0;

    fmpz_init_set_si(n, 138);

    flint_printf("Computing aliquot sequence\n");
    flint_printf("%3wu : ", 0);
    fmpz_print(n);
    flint_printf("\n");

    TIMEIT_ONCE_START

    while (!fmpz_is_one(n))
    {
        fmpz_sum_divisors_proper(n, n);
        num_iterations++;

        flint_printf("%3wu : ", num_iterations);
        fmpz_print(n);
        flint_printf("\n");
    }

    flint_printf("Sequence terminated after %wu iterations.\n", num_iterations);

    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    fmpz_clear(n);
    return 0;
}
