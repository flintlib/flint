/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

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
        printf("primegen N - print all primes <= N\n");
        printf("primegen -c N - generate the primes but just count them\n");
        return EXIT_FAILURE;
    }

    N = strtoul(argv[argc-1], NULL, 10);
    if (N == ULONG_MAX)
    {
        printf("N must be smaller than %lu\n", ULONG_MAX);
        return EXIT_FAILURE;
    }

    if (argc == 3)
    {
        ulong count = 0;
        n_primes_init(iter);
        while ((p = n_primes_next(iter)) <= N)
            count++;
        n_primes_clear(iter);
        printf("pi(%lu) = %lu\n", N, count);
    }
    else
    {
        n_primes_init(iter);
        while ((p = n_primes_next(iter)) <= N)
            printf("%lu\n", p);
        n_primes_clear(iter);
    }

    return EXIT_SUCCESS;
}
