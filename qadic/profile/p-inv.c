/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Benchmarks for the q-adic inversion routine.

    We consider the set-up with p = 17, N = 2^i, i = 0, ..., 19, 
    and invert A mod p^N, where 

        A = [a{0},...,a{d-1}], where a{i} = (3+i)^{3N}.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "qadic.h"

int
main(void)
{
    slong l, len = 20;
    slong runs[] = {
        1000, 1000, 1000, 1000, 1000, 
        1000, 100, 100, 100, 100, 
        100, 10, 10, 1, 1, 
        1, 1, 1, 1, 1
    };
    slong N[] = {
        1, 2, 4, 8, 16, 
        32, 64, 128, 256, 512, 
        1024, WORD(1) << 11, WORD(1) << 12, WORD(1) << 13, WORD(1) << 14, 
        WORD(1) << 15, WORD(1) << 16, WORD(1) << 17, WORD(1) << 18, WORD(1) << 19
    };
    slong T[20] = {0};

    flint_printf("Benchmark for q-adic inversion.\n");
    fflush(stdout);

for (l = 0; l < len; l++)
{
    FLINT_TEST_INIT(state);
    slong d = 97, i, n = N[l], r;
    clock_t c0, c1;
    long double cputime;

    fmpz_t p;
    qadic_ctx_t ctx;
    qadic_t a, z;

    

    fmpz_init_set_ui(p, 17);

    qadic_ctx_init_conway(ctx, p, d, n, n, "X", PADIC_VAL_UNIT);

    qadic_init(a);
    qadic_init(z);

    padic_poly_fit_length(a, d);
    _padic_poly_set_length(a, d);
    a->val = 0;

    for (i = 0; i < d; i++)
    {
        fmpz_t f, pow;

        fmpz_init(f);
        fmpz_init(pow);
        fmpz_set_ui(f, 3 + i);
        fmpz_pow_ui(pow, p, n);
        fmpz_pow_ui(a->coeffs + i, f, 3 * n);
        fmpz_mod(a->coeffs + i, a->coeffs + i, pow);
        fmpz_clear(f);
        fmpz_clear(pow);
    }
    _padic_poly_normalise(a);

    c0 = clock();
    for (r = runs[l]; (r); r--)
    {
        qadic_inv(z, a, ctx);
        qadic_zero(z);
    }
    c1 = clock();

    cputime = (long double) (c1 - c0) / (long double) CLOCKS_PER_SEC;

    T[l] = (slong) (cputime * (1000000000 / runs[l]));

    flint_printf("%2ld, %4XYXYXYXY, %8ld, %wd\n", 
        l, cputime, runs[l], T[l]);

    qadic_clear(a);
    qadic_clear(z);

    fmpz_clear(p);
    qadic_ctx_clear(ctx);
    flint_randclear(state);
}

    flint_printf("Output as a list:\n");
    for (l = 0; l < len; l++)
        flint_printf("%wd, ", T[l]);
    flint_printf("\n");
}

