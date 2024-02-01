/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    FLINT program for demonstrating the Integer Partition function.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"
#include "profiler.h"

int
main(int argc, char * argv[])
{
    fmpz_t x;
    ulong n;
    slong i;
    int quiet = 0;
    slong num_threads = 1;

    if (argc < 2)
    {
        flint_printf("usage: partitions n [-quiet]\n");
        return 1;
    }

    for (i = 2; i < argc; i++)
    {
        if (!strcmp(argv[i], "-quiet"))
            quiet = 1;
        else if (!strcmp(argv[i], "-threads"))
            num_threads = atol(argv[i+1]);
    }

    flint_set_num_threads(num_threads);

    flint_sscanf(argv[1], "%wu", &n);

    flint_printf("p(%wu) = \n", n);

    fmpz_init(x);
    TIMEIT_ONCE_START
    arith_number_of_partitions(x, n);
    TIMEIT_ONCE_STOP
    if (!quiet)
    {
        fmpz_print(x);
        flint_printf("\n");
    }
    fmpz_clear(x);

    return 0;
}

