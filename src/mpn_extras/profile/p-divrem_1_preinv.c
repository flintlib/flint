/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "profiler.h"

int main(void)
{
    slong N = 1000;
    slong maxn = 20;
    double t0, t1, FLINT_SET_BUT_UNUSED(tt);

    ulong d, dinv, norm;
    mp_ptr A, Q, R;
    slong n, i;

    A = flint_malloc(sizeof(mp_limb_t) * maxn * N);
    Q = flint_malloc(sizeof(mp_limb_t) * maxn * N);
    R = flint_malloc(sizeof(mp_limb_t) * N);

    mpn_random(A, maxn * N);

    for (norm = 0; norm <= 4; norm += 4)
    {
        flint_printf("n          mpn_divrem_1    preinv    speedup\n");

        d = n_nextprime(UWORD(1) << (FLINT_BITS - 1 - norm), 1);
        norm = flint_clz(d);
        dinv = n_preinvert_limb_prenorm(d << norm);
        flint_printf("\nd = %wu\n\n", d);

        for (n = 1; n <= maxn; n++)
        {
            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = mpn_divrem_1(Q + n * i, 0, A + n * i, n, d);
            TIMEIT_STOP_VALUES(tt, t0);
            TIMEIT_START;
            for (i = 0; i < N; i++)
            {
                R[i] = flint_mpn_divrem_1_preinv(Q + n * i, A + n * i, n, d, dinv, norm);
            }
            TIMEIT_STOP_VALUES(tt, t1);
            t0 /= N; t1 /= N;
            flint_printf("%10wd    %8g    %8g    %.3f\n", n, t0, t1, t0  / t1);
        }

        if (norm == 0)
        {
            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = mpn_divrem_1(Q + 2 * i, 0, A + 2 * i, 2, d);
            TIMEIT_STOP_VALUES(tt, t0);
            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = flint_mpn_divrem_2_1_preinv_norm(Q + 2 * i, A + 2 * i, d, dinv);
            TIMEIT_STOP_VALUES(tt, t1);
            t0 /= N; t1 /= N;
            flint_printf("    2_norm    %8g    %8g    %.3f\n", t0, t1, t0  / t1);

            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = mpn_divrem_1(Q + 3 * i, 0, A + 3 * i, 3, d);
            TIMEIT_STOP_VALUES(tt, t0);
            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = flint_mpn_divrem_3_1_preinv_norm(Q + 3 * i, A + 3 * i, d, dinv);
            TIMEIT_STOP_VALUES(tt, t1);
            t0 /= N; t1 /= N;
            flint_printf("    3_norm    %8g    %8g    %.3f\n", t0, t1, t0  / t1);
        }
        else
        {
            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = mpn_divrem_1(Q + i, 0, A + i, 1, d);
            TIMEIT_STOP_VALUES(tt, t0);
            TIMEIT_START;
            for (i = 0; i < N; i++)
            {
                R[i] = n_divrem_preinv_unnorm(Q + i, A[i], d, dinv, norm);
            }
            TIMEIT_STOP_VALUES(tt, t1);
            t0 /= N; t1 /= N;
            flint_printf("  1_unnorm    %8g    %8g    %.3f\n", t0, t1, t0  / t1);

            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = mpn_divrem_1(Q + 2 * i, 0, A + 2 * i, 2, d);
            TIMEIT_STOP_VALUES(tt, t0);
            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = flint_mpn_divrem_2_1_preinv_unnorm(Q + 2 * i, A + 2 * i, d, dinv, norm);
            TIMEIT_STOP_VALUES(tt, t1);
            t0 /= N; t1 /= N;
            flint_printf("  2_unnorm    %8g    %8g    %.3f\n", t0, t1, t0  / t1);

            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = mpn_divrem_1(Q + 3 * i, 0, A + 3 * i, 3, d);
            TIMEIT_STOP_VALUES(tt, t0);
            TIMEIT_START;
            for (i = 0; i < N; i++)
                R[i] = flint_mpn_divrem_3_1_preinv_unnorm(Q + 3 * i, A + 3 * i, d, dinv, norm);
            TIMEIT_STOP_VALUES(tt, t1);
            t0 /= N; t1 /= N;
            flint_printf("  3_unnorm    %8g    %8g    %.3f\n", t0, t1, t0  / t1);
        }

        flint_printf("\n");
    }

    flint_free(A);
    flint_free(Q);
    flint_free(R);

    flint_cleanup_master();
    return 0;
}

