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

            for (int check_reps = 0; check_reps < 2+50/(1+L); check_reps++)
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

            for (int check_reps = 0; check_reps < 2+90/(1+L); check_reps++)
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

/*
    Check that sd_fft_trunc and sd_ifft_trunc give correct results for inputs
    at the edges of their documented ranges: (-3n, 3n) for the forward
    transform and (-2n, 2n) for the inverse transform (the pointwise products
    fed to the inverse transform are in (-3/2 n, 3/2 n)). Random inputs in
    [0, n) almost never exercise the worst-case bounds of the butterflies, so
    this uses extreme values and compares with the transform of the same data
    reduced to [0, n); both transforms are linear mod n.
*/
static double sd_fft_test_extreme(flint_rand_t state, double lim, ulong i, int mode)
{
    if (mode == 0)          /* uniform in [-lim, lim] */
        return (double) n_randint(state, (ulong) (2*lim) + 1) - lim;
    else if (mode == 1)     /* +-lim, random signs */
        return n_randint(state, 2) ? lim : -lim;
    else                    /* +-lim, sign alternating in blocks of 1, 2, 4 or 8 */
        return ((i >> (mode - 2)) & 1) ? lim : -lim;
}

void test_sd_fft_trunc_extreme(sd_fft_ctx_t Q, ulong minL, ulong maxL, ulong reps, flint_rand_t state)
{
    for (ulong L = minL; L <= maxL; L++)
    {
        ulong i, Xn = n_pow2(L);
        double* a = (double*) flint_aligned_alloc(32, FLINT_MAX(32, Xn*sizeof(double)));
        double* b = (double*) flint_aligned_alloc(32, FLINT_MAX(32, Xn*sizeof(double)));

        sd_fft_ctx_fit_depth(Q, L);

        for (ulong rep = 0; rep < reps; rep++)
        {
            int mode = n_randint(state, 6);
            int inverse = n_randint(state, 2);
            /* keep the data at most 3n - 1 resp. 2n - 1 in absolute value */
            double lim = inverse ? 2*Q->p - 1 : 3*Q->p - 1;
            ulong itrunc = n_randint(state, 2) ? Xn : 1 + n_randint(state, Xn);
            ulong otrunc = n_randint(state, 2) ? Xn : 1 + n_randint(state, Xn);

            if (inverse)
                itrunc = otrunc;

            for (i = 0; i < Xn; i++)
            {
                a[i] = sd_fft_test_extreme(state, lim, i, mode);
                /* exact: |a[i]| < 2^53 */
                b[i] = a[i] - Q->p*floor(a[i]*(1.0/Q->p));
                while (b[i] < 0) b[i] += Q->p;
                while (b[i] >= Q->p) b[i] -= Q->p;
            }

            if (inverse)
            {
                sd_ifft_trunc(Q, a, L, otrunc);
                sd_ifft_trunc(Q, b, L, otrunc);
            }
            else
            {
                sd_fft_trunc(Q, a, L, itrunc, otrunc);
                sd_fft_trunc(Q, b, L, itrunc, otrunc);
            }

            for (i = 0; i < otrunc; i++)
            {
                if (a[i] != floor(a[i]) || !vec1d_same_mod(a[i], b[i], Q->p, Q->pinv))
                {
                    flint_printf("FAIL: %s with extreme inputs\n"
                                 "p: %wu\ndepth: %wu\nitrunc: %wu\notrunc: %wu\n"
                                 "mode: %d\nindex: %wu\ngot: %.1f\nexpected: %.1f (mod p)\n",
                                 inverse ? "ifft" : "fft", Q->mod.n, L, itrunc, otrunc,
                                 mode, i, a[i], b[i]);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        flint_aligned_free(a);
        flint_aligned_free(b);
    }
}

TEST_FUNCTION_START(sd_fft, state)
{
    {
        sd_fft_ctx_t Q;
        sd_fft_ctx_init_prime(Q, UWORD(0x0003f00000000001));
        test_sd_fft_trunc(Q, 0, 19, 20, state);
        test_sd_fft_trunc_extreme(Q, 0, 14, 40, state);
        sd_fft_ctx_clear(Q);
    }

    {
        /* two more of the default mpn_ctx primes */
        sd_fft_ctx_t Q;
        sd_fft_ctx_init_prime(Q, UWORD(1086317488242689));
        test_sd_fft_trunc_extreme(Q, 0, 14, 40, state);
        sd_fft_ctx_clear(Q);
        sd_fft_ctx_init_prime(Q, UWORD(659706976665601));
        test_sd_fft_trunc_extreme(Q, 0, 14, 40, state);
        sd_fft_ctx_clear(Q);
    }

    {
        sd_fft_ctx_t Q;
        sd_fft_ctx_init_prime(Q, UWORD(257));
        test_sd_fft_trunc(Q, 0, 8, 5, state);
        test_sd_fft_trunc_extreme(Q, 0, 8, 40, state);
        sd_fft_ctx_clear(Q);
    }

    TEST_FUNCTION_END(state);
}
