/*
    Copyright (C) 2026 Alexey Orlov

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Benchmarks for the q-adic square root routine with q=2^d.
 */

#include <time.h>

#include "fmpz.h"
#include "qadic.h"

#define BENCH_COUNT 3

int
main(void)
{
    const int test_deg[BENCH_COUNT] = {
        300, 500, 1000,
    };
    const int test_run[BENCH_COUNT] = {
        500, 500, 500
    };

    fmpz_t p;
    fmpz_init_set_ui(p, 2);

    int i, ci, r;

    for (i = 0; i < BENCH_COUNT; i++)
    {
        const int d = test_deg[i];
        const int R = test_run[i];

        qadic_ctx_t ctx;
        struct qadic2_sqrt_precomp *precomp;
        qadic_t a, b, z;

        clock_t c0, c1;

        qadic_ctx_init(ctx, p, d, 0, 1, "x", PADIC_SERIES);

        qadic_init(a);
        qadic_init(b);
        qadic_init(z);

        padic_poly_fit_length(a, d);
        _padic_poly_set_length(a, d);
        a->val = 0;

        fmpz_set_ui(a->coeffs + 0, 1);
        for (ci = 1; ci < d; ci++)
            fmpz_set_ui(a->coeffs + ci, 8);

        _padic_poly_normalise(a);
        qadic_mul(b, a, a, ctx);

        c0 = clock();
        for (r = 0; r < R; r++)
        {
            qadic_sqrt(z, b, ctx);
            qadic_zero(z);
        }
        c1 = clock();

        flint_printf("sqrt [d = %d, R = %d]: %.4Lf s\n", d, R, (long double)(c1 - c0) / (long double)CLOCKS_PER_SEC);

        precomp = _qadic_char2_sqrt_precomp_init(ctx);
        c0 = clock();
        for (r = 0; r < R; r++)
        {
            _qadic_char2_sqrt_with_precomp(z, b, ctx, precomp);
            qadic_zero(z);
        }
        c1 = clock();

        flint_printf("sqrt (precomp) [d = %d, R = %d]: %.4Lf s\n", d, R, (long double)(c1 - c0) / (long double)CLOCKS_PER_SEC);

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(z);

        _qadic_char2_sqrt_precomp_clear(precomp);
        qadic_ctx_clear(ctx);
    }

    fmpz_clear(p);

    return 0;
}

