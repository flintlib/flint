/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Benchmarks for the p-adic exponential method, rectangular.

    We consider the set-up with p = 17, N = 2^i, i = 0, ..., 19, 
    and compute the logarithm of e = 17 a + 1, a = 3^{3 N} mod p^N.
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
    long l, len = 20;
    long runs[] = {
        10000000, 1000000, 1000000, 100000, 100000, 
        10000, 10000, 10000, 1000, 1000, 
        100, 100, 10, 1, 1, 
        1, 1, 1, 1, 1
    };
    long N[] = {
        1, 2, 4, 8, 16, 
        32, 64, 128, 256, 512, 
        1024, WORD(1) << 11, WORD(1) << 12, WORD(1) << 13, WORD(1) << 14, 
        WORD(1) << 15, WORD(1) << 16, WORD(1) << 17, WORD(1) << 18, WORD(1) << 19
    };
    long T[20] = {0};

    flint_printf("Benchmark for p-adic logarithm (rectangular).\n");
    fflush(stdout);

for (l = 0; l < FLINT_MIN(16, len); l++)
{
    FLINT_TEST_INIT(state);
    long n = N[l], r;
    clock_t c0, c1;
    long double cputime;

    fmpz_t p;
    padic_ctx_t ctx;
    padic_t e, z;

    

    fmpz_init_set_ui(p, 17);

    padic_ctx_init(ctx, p, n, n, PADIC_VAL_UNIT);

    padic_init(e);
    padic_init(z);

    {
        fmpz_t f = {WORD(3)}, pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, p, n);
        fmpz_pow_ui(padic_unit(e), f, 3 * n);
        fmpz_mul_ui(padic_unit(e), padic_unit(e), 17);
        fmpz_add_ui(padic_unit(e), padic_unit(e), 1);
        fmpz_mod(padic_unit(e), padic_unit(e), pow);
        fmpz_clear(pow);
    }

    c0 = clock();
    for (r = runs[l]; (r); r--)
    {
        padic_log_rectangular(z, e, ctx);
        padic_zero(z);
    }
    c1 = clock();

    padic_log_rectangular(z, e, ctx);
    padic_exp_rectangular(z, z, ctx);
    if (!padic_equal(e, z))
    {
        flint_printf("FAIL:\n");
        flint_printf("e = "), padic_print(e, ctx), flint_printf("\n");
        flint_printf("z = "), padic_print(z, ctx), flint_printf("\n");
        flint_printf("p = %wd\n", *p), flint_printf("\n");
        flint_printf("N = %wd\n", n), flint_printf("\n");
        flint_abort();
    }

    cputime = (long double) (c1 - c0) / (long double) CLOCKS_PER_SEC;

    T[l] = (slong) (cputime * (1000000000 / runs[l]));

    flint_printf("%2ld, %4XYXYXYXY, %8ld, %wd\n", 
        l, cputime, runs[l], T[l]);

    padic_clear(e);
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

