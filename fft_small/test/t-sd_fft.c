/* 
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"
#include "machine_vectors.h"

vec1d vec1d_eval_poly_mod(const vec1d* a, ulong an, const vec1d b, const vec1d n, const vec1d ninv)
{
    vec1d x = a[--an];
    while (an > 0)
        x = vec1d_add(a[--an], vec1d_mulmod2(x, b, n, ninv));
    return vec1d_reduce_to_pm1n(x, n, ninv);
}


void flint_print_d_fixed(double x, ulong l)
{
    ulong i;
    TMP_INIT;
    TMP_START;
    char* s = TMP_ARRAY_ALLOC(l + 1, char);

    for (i = 0; i < l; i++)
        s[i] = ' ';
    s[l] = 0;

    ulong y = fabs(rint(x));
    while (l > 0)
    {
        s[--l] = '0' + (y%10);
        y = y/10;
        if (y == 0)
            break;
    }

    printf("%s", s);

    TMP_END;
}


void test_v2_fft(sd_fft_ctx_t Q, ulong minL, ulong maxL, ulong ireps, flint_rand_t state)
{
    ulong irepmul = 10;
    minL = n_max(minL, LG_BLK_SZ);

    ulong dgs1 = n_sizeinbase(maxL, 10);
    ulong dgs2 = n_sizeinbase(ireps + irepmul*maxL, 10);

    flint_printf(" fft ");
    flint_print_d_fixed(0, dgs1);
    flint_printf(".");
    flint_print_d_fixed(0, dgs2);
    flint_printf("/");
    flint_print_d_fixed(maxL, dgs1);
    flint_printf(".");
    flint_print_d_fixed(0, dgs2);
    fflush(stdout);

    for (ulong L = minL; L <= maxL; L++)
    {
        ulong i;
        ulong Xn = n_pow2(L);
        double* X = FLINT_ARRAY_ALLOC(Xn, double);
        double* data =  (double*) flint_aligned_alloc(32,
                                      sd_fft_ctx_data_size(L)*sizeof(double));

        ulong nreps = ireps + irepmul*L;
        for (ulong rep = 0; rep < nreps; rep++)
        {
            for (ulong ii = 0; ii < 2*(dgs1 + dgs2) + 8; ii++)
                flint_printf("%c", '\b');

            flint_printf(" fft ");
            flint_print_d_fixed(L, dgs1);
            flint_printf(".");
            flint_print_d_fixed(1+rep, dgs2);
            flint_printf("/");
            flint_print_d_fixed(maxL, dgs1);
            flint_printf(".");
            flint_print_d_fixed(nreps, dgs2);
            fflush(stdout);

            // randomize input data
            for (i = 0; i < Xn; i++)
                X[i] = i + 1;

            // output of fft_trunc is supposed to be eval_poly
            ulong itrunc = n_round_up(1 + n_randint(state, Xn), Q->blk_sz);
            ulong otrunc = n_round_up(1 + n_randint(state, Xn) , Q->blk_sz);

            for (i = 0; i < itrunc; i++)
                sd_fft_ctx_set_index(data, i, X[i]);

            sd_fft_ctx_fft_trunc(Q, data, L, itrunc, otrunc);

            for (int check_reps = 0; check_reps < 3; check_reps++)
            {
                i = n_randint(state, otrunc);
                double point = sd_fft_ctx_w(Q, i);
                double y = vec1d_eval_poly_mod(X, itrunc, point, Q->p, Q->pinv);
                if (!vec1d_same_mod(y, sd_fft_ctx_get_fft_index(data, i), Q->p, Q->pinv))
                {
                    flint_printf("FAIL: fft error at index %wu\nitrunc: %wu\n"
                                           "otrunc: %wu\n", i, itrunc, otrunc);
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* output of ifft_trunc is supposed to be 2^L*input */
            ulong trunc = n_round_up(1 + n_randint(state, Xn), Q->blk_sz);
            for (i = 0; i < trunc; i++)
                sd_fft_ctx_set_index( data, i, X[i]);

            sd_fft_ctx_fft_trunc(Q, data, L, trunc, trunc);
            sd_fft_ctx_ifft_trunc(Q, data, L, trunc);

            for (int check_reps = 0; check_reps < 3; check_reps++)
            {
                i = n_randint(state, trunc);
                double m = vec1d_reduce_0n_to_pmhn(nmod_pow_ui(2, L, Q->mod), Q->p);
                double y = vec1d_mulmod2(X[i], m, Q->p, Q->pinv);
                if (!vec1d_same_mod(y, sd_fft_ctx_get_index(data, i), Q->p, Q->pinv))
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

    for (ulong ii = 0; ii < 2*(dgs1 + dgs2) + 8; ii++)
        flint_printf("%c", '\b');
    fflush(stdout);
}


int main(void)
{
    FLINT_TEST_INIT(state);

    flint_printf("sd_fft....");
    fflush(stdout);

    {
        sd_fft_ctx_t Q;
        sd_fft_ctx_init_prime(Q, UWORD(0x0003f00000000001));
        test_v2_fft(Q, 10, 19, 20, state);
        sd_fft_ctx_clear(Q);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
