/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "flint.h"
#include "ulong_extras.h"

int main(int argc, char* argv[])
{
    n_primes_t iter;
    mp_limb_t p, N;

    if (argc < 2)
    {
        flint_printf("primegen N - print all primes <= N\n");
        flint_printf("primegen -c N - generate the primes but just count them\n");
        return EXIT_FAILURE;
    }

    N = strtoul(argv[argc-1], NULL, 10);
    if (N == UWORD_MAX)
    {
        flint_printf("N must be smaller than %wu\n", UWORD_MAX);
        return EXIT_FAILURE;
    }

    if (argc == 3)
    {
        ulong count = 0;
        n_primes_init(iter);
        while ((p = n_primes_next(iter)) <= N)
            count++;
        n_primes_clear(iter);
        flint_printf("pi(%wu) = %wu\n", N, count);
    }
    else
    {
        n_primes_init(iter);
        while ((p = n_primes_next(iter)) <= N)
            flint_printf("%wu\n", p);
        n_primes_clear(iter);
    }

    return EXIT_SUCCESS;
}
