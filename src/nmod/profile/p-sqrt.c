/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "ulong_extras.h"
#include "nmod.h"
#include "profiler.h"

#define N 1000

/* random primes of the given bit size, either arbitrary or with
   n - 1 = q 2^s, q odd, for given s */
static void
random_primes(nn_ptr P, nmod_t * M, nn_ptr A, flint_rand_t state, int bits, int s, int squares)
{
    slong i;
    ulong p, c;

    for (i = 0; i < N; i++)
    {
        if (s == 0)
        {
            p = n_randprime(state, bits, 0);
        }
        else
        {
            do
            {
                c = (n_randlimb(state) >> (FLINT_BITS - (bits - s))) | 1;
                c |= UWORD(1) << (bits - s - 1);
                p = (c << s) + 1;
            }
            while (!n_is_prime(p));
        }

        P[i] = p;
        nmod_init(&M[i], p);
        A[i] = n_randlimb(state) % p;
        if (squares)
            A[i] = nmod_mul(A[i], A[i], M[i]);
    }
}

int main(void)
{
    flint_rand_t state;
    ulong P[N], A[N], r, sum = 0;
    nmod_t M[N];
    int bits, s, squares;
    slong i;
    double t1, FLINT_SET_BUT_UNUSED(tcpu);

    int bits_arr[] = { 8, 16, 24, 30, 32, 40, 50, 60, 62, 63, 64, 0 };
    int s_arr[] = { 0, 3, 8, 16, 24, 32, 48, -1 };
    int bi, si;

    flint_rand_init(state);

    flint_printf("Time (ns) per call; s = 0 means random primes\n");
    flint_printf("%5s %3s %12s %12s %12s\n", "bits", "s", "n_sqrtmod", "nmod_sqrt", "nmod_sqrt sq");

    for (bi = 0; (bits = bits_arr[bi]) != 0; bi++)
    {
        if (bits > FLINT_BITS)
            continue;

        for (si = 0; (s = s_arr[si]) >= 0; si++)
        {
            if (s > bits - 4)
                continue;

            flint_printf("%5d %3d", bits, s);

            for (squares = 0; squares <= 2; squares++)
            {
                random_primes(P, M, A, state, bits, s, squares == 2);

                TIMEIT_START;
                if (squares == 0)
                {
                    for (i = 0; i < N; i++)
                        sum += n_sqrtmod(A[i], P[i]);
                }
                else
                {
                    for (i = 0; i < N; i++)
                    {
                        nmod_sqrt(&r, A[i], M[i]);
                        sum += r;
                    }
                }
                TIMEIT_STOP_VALUES(tcpu, t1);

                flint_printf(" %12.1f", t1 * 1e9 / N);
            }

            flint_printf("\n");
        }
    }

    if (sum == 12345)
        flint_printf("\n");

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
