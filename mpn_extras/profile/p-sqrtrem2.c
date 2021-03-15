/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include "profiler.h"
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

double clocktime()
{
    return (double) clock() / CLOCKS_PER_SEC;
}

int main()
{
    flint_rand_t state;
    flint_randinit(state);
    slong i;

    slong N = 1000000;
    slong bits;

    printf("bits   mpn_sqrtrem   flint_mpn_sqrtrem2   speedup\n");

    for (bits = 128; bits >= 65; bits--)
    {
        mp_ptr X, S1, S2, R1, R2;
        int *L1, *L2;
        double t1, t2;

        X = flint_malloc(2 * N * sizeof(mp_limb_t));
        S1 = flint_malloc(N * sizeof(mp_limb_t));
        S2 = flint_malloc(N * sizeof(mp_limb_t));
        R1 = flint_malloc(2 * N * sizeof(mp_limb_t));
        R2 = flint_malloc(2 * N * sizeof(mp_limb_t));
        L1 = flint_malloc(N * sizeof(int));
        L2 = flint_malloc(N * sizeof(int));

        for (i = 0; i < N; i++)
        {
#if 0
            X[2 * i] = n_randtest(state);
            X[2 * i + 1] = n_randtest(state);
#else
            X[2 * i] = n_randlimb(state);
            X[2 * i + 1] = n_randlimb(state);
#endif

            X[2 * i + 1] |= (UWORD(1) << 63);
            X[2 * i + 1] >>= (128 - bits);

#if 0
            /* test perfect squares, or +/- 1 */
            if (n_randint(state, 2))
            {
                mpn_sqrtrem(S1, NULL, X + 2 * i, 2);
                umul_ppmm(X[2 * i + 1], X[2 * i], S1[0], S1[0]);
                if (n_randint(state, 2))
                    add_ssaaaa(X[2 * i + 1], X[2 * i], X[2 * i + 1], X[2 * i], 0, 1);
                if (n_randint(state, 2))
                    sub_ddmmss(X[2 * i + 1], X[2 * i], X[2 * i + 1], X[2 * i], 0, 1);
            }
#endif
        }

        t1 = clocktime();
        for (i = 0; i < N; i++)
        {
            L1[i] = mpn_sqrtrem(S1 + i, R1 + 2 * i, X + 2 * i, 2);
        }
        t1 = (clocktime() - t1) / N;

        t2 = clocktime();
        for (i = 0; i < N; i++)
        {
            L2[i] = flint_mpn_sqrtrem2(S2 + i, R2 + 2 * i, X + 2 * i);
        }
        t2 = (clocktime() - t2) / N;

        for (i = 0; i < N; i++)
        {
            if (S1[i] != S2[i] || L1[i] != L2[i] || (L1[i] > 0 && R1[2 * i] != R2[2 * i]) || (L1[i] == 2 && R1[2 * i + 1] != R2[2 * i + 1]))
            {
                printf("FAIL X = %lu %lu  S = %lu %lu  Rhi = %lu %lu  Rlo = %lu %lu\n",
                    X[2 * i], X[2 * i + 1], S1[i], S2[i], R1[2 * i + 1], R2[2 * i + 1], R1[2 * i], R2[2 * i]);
                abort();
            }
        }

        printf("%3ld   %g   %g   %f\n", bits, t1, t2, t1 / t2);

        free(X);
        free(S1);
        free(S2);
        free(R1);
        free(R2);
        free(L1);
        free(L2);
    }

    return 0;
}
