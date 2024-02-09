/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "profiler.h"

#define MAXN 16
#define MAXN2 256

int main(void)
{
    flint_rand_t state;
    mp_limb_t x[MAXN2], y[MAXN2], r[2 * MAXN2], s[2 * MAXN2];
    slong i, n, m;
    double t1, t2, __attribute__((unused)) __;

    flint_randinit(state);

    {
        for (i = 0; i < MAXN; i++)
        {
            x[i] = n_randtest(state);
            y[i] = n_randtest(state);
        }

        flint_printf("mpn_mul_n vs flint_mpn_mul_n\n\n");

        for (n = 1; n <= MAXN2; n = (n <= MAXN ? n + 1 : n * 1.2))
        {
            flint_printf("n = %wd    ", n);

            TIMEIT_START
#if 0
            __gmpn_mul_basecase(r, x, n, y, n);
#else
            mpn_mul_n(r, x, y, n);
#endif
            TIMEIT_STOP_VALUES(__, t1)
            TIMEIT_START
            flint_mpn_mul_n(s, x, y, n);
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%g    %g   %.2f\n", t1, t2, t1 / t2);

            if (mpn_cmp(r, s, 2 * n) != 0)
                flint_abort();
        }

        flint_printf("mpn_mul vs flint_mpn_mul\n\n");

        for (n = 1; n <= MAXN2; n = (n <= MAXN ? n + 1 : n * 1.2))
        {
            flint_printf("n = %wd    ", n);

            for (m = 1; m <= n; m = (m <= MAXN ? m + 1 : m * 1.2))
            {
                TIMEIT_START
                mpn_mul(r, x, n, y, m);
                TIMEIT_STOP_VALUES(__, t1)
                TIMEIT_START
                flint_mpn_mul(s, x, n, y, m);
                TIMEIT_STOP_VALUES(__, t2)

                flint_printf("%.2f ", t1 / t2);

                if (mpn_cmp(r, s, n + m) != 0)
                    flint_abort();
            }

            flint_printf("\n");
        }
    }

    flint_randclear(state);
    flint_cleanup_master();
    return 0;
}
