/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    For each depth and truncation, times sd_fft_trunc, sd_ifft_trunc and their
    transposes sd_fft_trunc_tr and sd_ifft_trunc_tr, in "gflops" counting
    5*len*log2(len) operations per transform of length len.

    A transposed transform performs exactly the same additions and
    multiplications as the transform it transposes, in the reverse order, so
    the two columns of each pair should be close; the "ratio" columns give
    the transposed time over the direct time.

    Usage:
        build/fft_small/profile/p-sd_fft [mindepth [maxdepth]]
*/

#include <stdlib.h>
#include "ulong_extras.h"
#include "fft_small.h"
#include "profiler.h"

static void print_d_fixed_dot(double x, ulong l, ulong r)
{
    ulong i;
    TMP_INIT;
    TMP_START;
    char* s = TMP_ARRAY_ALLOC(l + r + 2, char);

    for (i = 0; i < l + r + 1; i++)
        s[i] = ' ';
    s[l + r + 1] = 0;

    s[l] = '.';
    ulong y = fabs(rint(x*n_pow(10, r)));
    while (r > 0)
    {
        --r;
        s[l+1+r] = '0' + (y%10);
        y = y/10;
    }
    while (l > 0)
    {
        --l;
        s[l] = '0' + (y%10);
        y = y/10;
        if (y == 0)
            break;
    }

    printf("%s", s);

    TMP_END;
}

static void record_time(int first, double len, double t_f, double t_ft,
                                               double t_i, double t_it)
{
    double depth = log2(len);

    t_f  *= 1e6;
    t_ft *= 1e6;
    t_i  *= 1e6;
    t_it *= 1e6;

    flint_printf(first ? "depth" : "     ");
    print_d_fixed_dot(depth, 3, 1);
    flint_printf(":");
    print_d_fixed_dot(5*len*depth/t_f,  4, 2);
    print_d_fixed_dot(5*len*depth/t_ft, 4, 2);
    print_d_fixed_dot(t_ft/t_f, 3, 2);
    flint_printf("   |");
    print_d_fixed_dot(5*len*depth/t_i,  4, 2);
    print_d_fixed_dot(5*len*depth/t_it, 4, 2);
    print_d_fixed_dot(t_it/t_i, 3, 2);
    flint_printf("\n");
}

#define EXP_FAC 16

static void profile_sd_fft(sd_fft_ctx_t Q, ulong minL, ulong maxL)
{
    timeit_t timer;

    minL = n_max(minL, 10);
    maxL = n_max(maxL, minL + 1);

    timeit_start(timer);
    sd_fft_ctx_fit_depth(Q, maxL);
    timeit_stop(timer);
    flint_printf("depth %wu setup: %wd ms\n", maxL, timer->wall);

    flint_printf("             fft   fft_tr ratio   |    ifft  ifft_tr ratio\n");

    minL = n_max(minL, LG_BLK_SZ + 1);
    for (ulong L = minL+1; L <= maxL; L++)
    {
        ulong  sz = sd_fft_ctx_data_size(L)*sizeof(double);
        double* data = flint_aligned_alloc(4096, n_round_up(sz, 4096));

        /* 1/2*2^L < otrunc <= 2^L */
        ulong otrunc = n_pow2(L-1);
        while (otrunc < n_pow2(L))
        {
            double t_f, t_ft, t_i, t_it;

            otrunc = n_round_up(1+1+(EXP_FAC+1)*otrunc/EXP_FAC, BLK_SZ);
            otrunc = n_min(otrunc, n_pow2(L));
            ulong itrunc = n_round_up(otrunc/2, BLK_SZ);
            for (ulong i = 0; i < n_pow2(L); i++)
                data[i] = 0;
            ulong nreps = 1 + 300000000/(otrunc*n_clog2(otrunc));

            timeit_start(timer);
            for (ulong i = 0; i < nreps; i++)
                sd_fft_trunc(Q, data, L, itrunc, otrunc);
            timeit_stop(timer);
            t_f = (double) timer->wall / nreps;

            timeit_start(timer);
            for (ulong i = 0; i < nreps; i++)
                sd_fft_trunc_tr(Q, data, L, itrunc, otrunc);
            timeit_stop(timer);
            t_ft = (double) timer->wall / nreps;

            timeit_start(timer);
            for (ulong i = 0; i < nreps; i++)
                sd_ifft_trunc(Q, data, L, otrunc);
            timeit_stop(timer);
            t_i = (double) timer->wall / nreps;

            timeit_start(timer);
            for (ulong i = 0; i < nreps; i++)
                sd_ifft_trunc_tr(Q, data, L, otrunc);
            timeit_stop(timer);
            t_it = (double) timer->wall / nreps;

            record_time(otrunc == n_pow2(L), otrunc, t_f, t_ft, t_i, t_it);
        }

        flint_free(data);
    }
}

int main(int argc, char** argv)
{
    sd_fft_ctx_t Q;
    ulong minL = (argc > 1) ? (ulong) atoi(argv[1]) : 9;
    ulong maxL = (argc > 2) ? (ulong) atoi(argv[2]) : 20;

    sd_fft_ctx_init_prime(Q, UWORD(0x0003f00000000001));
    flint_printf("--- fft/ifft and transposes, gflops ---\n");
    profile_sd_fft(Q, minL, maxL);
    sd_fft_ctx_clear(Q);

    return 0;
}
