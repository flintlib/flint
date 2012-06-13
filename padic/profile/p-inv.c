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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

/*
    Benchmarks for the p-adic inversion routine.

    We consider the set-up with p = 17, N = 2^i, i = 0, ..., 19, 
    and invert a = 3^{3 N} mod p^N.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "padic.h"

int
main(void)
{
    long l, len = 20;
    long runs[] = {
        10000000, 1000000, 1000000, 1000000, 1000000, 
        100000, 100000, 100000, 100000, 10000, 
        10000, 1000, 1000, 1000, 10, 
        10, 10, 10, 10, 1
    };
    long N[] = {
        1, 2, 4, 8, 16, 
        32, 64, 128, 256, 512, 
        1024, 1L << 11, 1L << 12, 1L << 13, 1L << 14, 
        1L << 15, 1L << 16, 1L << 17, 1L << 18, 1L << 19
    };
    long T[20] = {0};

    printf("Benchmark for p-adic inversion.\n");
    fflush(stdout);

for (l = 0; l < len; l++)
{
    flint_rand_t state;
    long n = N[l], r;
    clock_t c0, c1;
    long double cputime;

    fmpz_t p;
    padic_ctx_t ctx;
    padic_t a, z;

    flint_randinit(state);

    fmpz_init_set_ui(p, 17);

    padic_ctx_init(ctx, p, n, PADIC_VAL_UNIT);

    padic_init(a, ctx);
    padic_init(z, ctx);

    {
        fmpz_t f = {3L}, pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, p, n);
        fmpz_pow_ui(padic_unit(a), f, 3 * n);
        fmpz_mod(padic_unit(a), padic_unit(a), pow);
        fmpz_clear(pow);
    }

    c0 = clock();
    for (r = runs[l]; (r); r--)
    {
        padic_inv(z, a, ctx);
        padic_zero(z);
    }
    c1 = clock();

    cputime = (long double) (c1 - c0) / (long double) CLOCKS_PER_SEC;

    T[l] = (long) (cputime * (1000000000 / runs[l]));

    printf("%2ld, %4LG, %8ld, %ld\n", 
        l, cputime, runs[l], T[l]);

    padic_clear(a, ctx);
    padic_clear(z, ctx);

    fmpz_clear(p);
    padic_ctx_clear(ctx);
    flint_randclear(state);
}

    printf("Output as a list:\n");
    for (l = 0; l < len; l++)
        printf("%ld, ", T[l]);
    printf("\n");
}

