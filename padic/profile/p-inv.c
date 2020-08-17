/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "padic.h"

int
main(void)
{
    slong l, len = 20;
    slong runs[] = {
        10000000, 1000000, 1000000, 1000000, 1000000, 
        100000, 100000, 100000, 100000, 10000, 
        10000, 1000, 1000, 1000, 10, 
        10, 10, 10, 10, 1
    };
    slong N[] = {
        1, 2, 4, 8, 16, 
        32, 64, 128, 256, 512, 
        1024, WORD(1) << 11, WORD(1) << 12, WORD(1) << 13, WORD(1) << 14, 
        WORD(1) << 15, WORD(1) << 16, WORD(1) << 17, WORD(1) << 18, WORD(1) << 19
    };
    slong T[20] = {0};

    flint_printf("Benchmark for p-adic inversion.\n");
    fflush(stdout);

for (l = 0; l < len; l++)
{
    FLINT_TEST_INIT(state);
    slong n = N[l], r;
    clock_t c0, c1;
    long double cputime;

    fmpz_t p;
    padic_ctx_t ctx;
    padic_t a, z;

    

    fmpz_init_set_ui(p, 17);

    padic_ctx_init(ctx, p, n, n, PADIC_VAL_UNIT);

    padic_init(a);
    padic_init(z);

    {
        fmpz_t f = {WORD(3)}, pow;

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

    T[l] = (slong) (cputime * (1000000000 / runs[l]));

    flint_printf("%2ld, %4XYXYXYXY, %8ld, %wd\n", 
        l, cputime, runs[l], T[l]);

    padic_clear(a);
    padic_clear(z);

    fmpz_clear(p);
    padic_ctx_clear(ctx);
    flint_randclear(state);
}

    flint_printf("Output as a list:\n");
    for (l = 0; l < len; l++)
        flint_printf("%wd, ", T[l]);
    flint_printf("\n");
}

