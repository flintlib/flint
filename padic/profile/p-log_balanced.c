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
    Benchmarks for the p-adic exponential method, balanced.

    We consider the set-up with p = 17, N = 2^i, i = 0, ..., 19, 
    and compute the logarithm of e = 17 a + 1, a = 3^{3 N} mod p^N.
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
        10000000, 1000000, 100000, 100000, 100000, 
        10000, 10000, 10000, 1000, 1000, 
        1000, 100, 100, 10, 10, 
        1, 1, 1, 1, 1
    };
    long N[] = {
        1, 2, 4, 8, 16, 
        32, 64, 128, 256, 512, 
        1024, 1L << 11, 1L << 12, 1L << 13, 1L << 14, 
        1L << 15, 1L << 16, 1L << 17, 1L << 18, 1L << 19
    };
    long T[20] = {0};

    printf("Benchmark for p-adic logarithm (balanced).\n");
    fflush(stdout);

for (l = 0; l < len; l++)
{
    flint_rand_t state;
    long n = N[l], r;
    clock_t c0, c1;
    long double cputime;

    fmpz_t p;
    padic_ctx_t ctx;
    padic_t e, z;

    flint_randinit(state);

    fmpz_init_set_ui(p, 17);

    padic_ctx_init(ctx, p, n, PADIC_VAL_UNIT);

    padic_init(e, ctx);
    padic_init(z, ctx);

    {
        fmpz_t f = {3L}, pow;

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
        padic_log_balanced(z, e, ctx);
        padic_zero(z);
    }
    c1 = clock();

    padic_log_balanced(z, e, ctx);
    padic_exp_balanced(z, z, ctx);
    if (!_padic_equal(e, z))
    {
        printf("FAIL:\n");
        printf("e = "), padic_print(e, ctx), printf("\n");
        printf("z = "), padic_print(z, ctx), printf("\n");
        printf("p = %ld\n", *p), printf("\n");
        printf("N = %ld\n", n), printf("\n");
        abort();
    }

    cputime = (long double) (c1 - c0) / (long double) CLOCKS_PER_SEC;

    T[l] = (long) (cputime * (1000000000 / runs[l]));

    printf("%2ld, %4LG, %8ld, %ld\n", 
        l, cputime, runs[l], T[l]);

    padic_clear(e, ctx);
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

