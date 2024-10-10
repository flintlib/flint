/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "ulong_extras.h"
#include "fft_small.h"
#include "machine_vectors.h"

vec1d vec1d_eval_poly_mod(const vec1d* a, ulong an, const vec1d b, const vec1d n, const vec1d ninv)
{
    vec1d x = a[--an];
    while (an > 0)
        x = vec1d_add(a[--an], vec1d_mulmod(x, b, n, ninv));
    return vec1d_reduce_to_pm1n(x, n, ninv);
}

void test_sd_fft_trunc(sd_fft_ctx_t Q, ulong minL, ulong maxL, ulong ireps, flint_rand_t state)
{
    ulong irepmul = 10;

    for (ulong L = minL; L <= maxL; L++)
    {
        ulong i;
        ulong Xn = n_pow2(L);
        double* X = FLINT_ARRAY_ALLOC(Xn, double);
        double* data =  (double*) flint_aligned_alloc(32,
                                      FLINT_MAX(32, n_pow2(L)*sizeof(double)));

        ulong nreps = ireps + irepmul*L;
        for (ulong rep = 0; rep < nreps; rep++)
        {
            /* randomize input data */
            for (i = 0; i < Xn; i++)
                X[i] = n_randint(state, Q->mod.n);

            /* output of fft_trunc is supposed to be eval_poly */
            ulong itrunc = rep == 0 ? Xn : 1 + n_randint(state, Xn);
            ulong otrunc = rep == 0 ? Xn : 1 + n_randint(state, Xn);

            for (i = 0; i < itrunc; i++)
                data[i] = X[i];

            sd_fft_trunc(Q, data, L, itrunc, otrunc);

            for (ulong check_reps = 0; check_reps < 2+50/(1+L); check_reps++)
            {
                i = n_randint(state, otrunc);
                double point = sd_fft_ctx_w(Q, i);
                double y = vec1d_eval_poly_mod(X, itrunc, point, Q->p, Q->pinv);
                if (!vec1d_same_mod(y, data[sd_fft_ctx_trunc_index(L, i)], Q->p, Q->pinv))
                {
                    flint_printf("FAIL: fft error at index %wu\nitrunc: %wu\n"
                                           "otrunc: %wu\ndepth: %wu\n", i, itrunc, otrunc, L);
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* output of ifft_trunc is supposed to be 2^L*input */
            ulong trunc = rep == 0 ? Xn : 1 + n_randint(state, Xn);
            for (i = 0; i < trunc; i++)
                data[i] = X[i];

            sd_fft_trunc(Q, data, L, trunc, trunc);
            sd_ifft_trunc(Q, data, L, trunc);

            for (ulong check_reps = 0; check_reps < 2+90/(1+L); check_reps++)
            {
                i = n_randint(state, trunc);
                double m = vec1d_reduce_0n_to_pmhn(nmod_pow_ui(2, L, Q->mod), Q->p);
                double y = vec1d_mulmod(X[i], m, Q->p, Q->pinv);
                if (!vec1d_same_mod(y, data[i], Q->p, Q->pinv))
                {
                    flint_printf("FAIL: ifft error at index %wu\n"
                                      "trunc: %wu\ndepth: %wu\n", i, trunc, L);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        flint_aligned_free(data);
        flint_free(X);
    }
}

TEST_FUNCTION_START(sd_fft, state)
{
    {
        sd_fft_ctx_t Q;
        sd_fft_ctx_init_prime(Q, UWORD(0x0003f00000000001));
        test_sd_fft_trunc(Q, 0, 19, 20, state);
        sd_fft_ctx_clear(Q);
    }

    TEST_FUNCTION_END(state);
}
