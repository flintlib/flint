/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Tuning program for the per-prime-count weights npw[] used by
    _mulmod_2expp1_choose_m in fft_small/negacyclic.c.

    Usage:  p-mulmod_2expp1_tuning [big]

    Races every admissible (m, b, np) geometry of a set of
    Schoenhage-Strassen-style ring sizes, computes the cost per model
    unit np * m * (depth + 4 + b/16) for each, aggregates a median per
    prime count, and prints a paste-ready weight table normalized to
    sixteenths of the np = 3..4 baseline, followed by a check that the
    suggested table reproduces the measured winner at every raced size.
    Pass "big" to add two multi-million-bit sizes (slower).
*/

#include <stdio.h>
#include <string.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"
#include "fft_small.h"
#include "profiler.h"

#define TIMEIT(var, reps, code) \
    { timeit_t t_; double best_ = 1e30; \
      for (int r_ = 0; r_ < 4; r_++) { \
        timeit_start(t_); \
        for (slong i_ = 0; i_ < (reps); i_++) { code; } \
        timeit_stop(t_); \
        if ((double) t_->cpu / (reps) < best_) \
            best_ = (double) t_->cpu / (reps); } \
      var = 1000.0 * best_; }

typedef struct { slong N, m, b, np; double t, units; } sample_t;

/* replicate the chooser's admissibility: capacity for some np <= 8,
   and digits (b <= 50) or an instantiated conversion */
static int
probe(mpn_ctx_struct * R, slong b, slong depth)
{
    fmpz_t P;
    slong np;
    int ok = 0;

    fmpz_init_set_ui(P, R->ffts[0].p);
    for (np = 2; np <= 8 && !ok; np++)
    {
        fmpz_mul_ui(P, P, R->ffts[np - 1].p);
        if ((slong) fmpz_bits(P) - 1 >= 2 * b + depth + 1
            && (b <= 50
                || _mpn_ctx_to_ffts_func((ulong) np, (ulong) b) != NULL))
            ok = 1;
    }
    fmpz_clear(P);
    return ok;
}

int main(int argc, char ** argv)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    slong Ns[16] = {47104, 69632, 172032, 237568, 786432,
                    1179648, 1835008, 0};
    slong nN = 7;
    sample_t samp[256];
    slong nsamp = 0;
    double med[9];
    slong cnt[9] = {0};
    slong w[9] = {0};
    slong i, j, np;
    double base;

    if (argc > 1 && strcmp(argv[1], "big") == 0)
    {
        Ns[nN++] = 2752512;
        Ns[nN++] = 6029312;
    }

    flint_set_num_threads(1);
    flint_printf("racing admissible negacyclic geometries "
                 "(this takes a few minutes)\n\n");

    for (i = 0; i < nN; i++)
    {
        slong N = Ns[i], m;

        flint_printf("N %wd:\n", N);
        for (m = 16; m <= N / 20; m *= 2)
        {
            slong b, depth, h, reps, k;
            sd_fft_mpn_mulmod_2expp1_ctx_struct C;
            sd_fft_mpn_mulmod_2expp1_scratch_struct S;
            nn_ptr x, y, z;
            double * F;
            double t;

            if (N % m != 0)
                continue;
            b = N / m;
            if (b < 20 || b > 200)
                continue;
            depth = FLINT_BIT_COUNT(m) - 1;
            if (!probe(R, b, depth))
                continue;

            h = N / 64;
            sd_fft_mpn_mulmod_2expp1_ctx_init(&C, R, N, m);
            sd_fft_mpn_mulmod_2expp1_scratch_init(&S, &C);
            x = flint_malloc((h + 1) * sizeof(ulong));
            y = flint_malloc((h + 1) * sizeof(ulong));
            z = flint_malloc((h + 1) * sizeof(ulong));
            for (k = 0; k <= h; k++)
            {
                x[k] = 7 * k + 1;
                y[k] = 3 * k + 11;
            }
            x[h] &= 1;
            y[h] &= 1;
            F = flint_aligned_alloc(4096,
                ((sd_fft_mpn_mulmod_2expp1_transformed_size(&C)
                  * sizeof(double) + 4095) / 4096) * 4096);
            sd_fft_mpn_mulmod_2expp1_transform(&C, F, y, &S);
            reps = 30000000 / (h * (slong) FLINT_BIT_COUNT(h)) + 20;
            TIMEIT(t, reps,
                   sd_fft_mpn_mulmod_2expp1_mul_cached(&C, z, x, F, &S))

            samp[nsamp].N = N;
            samp[nsamp].m = C.m;
            samp[nsamp].b = C.b;
            samp[nsamp].np = C.np;
            samp[nsamp].t = t;
            samp[nsamp].units = (double) C.np * (double) C.m
                * ((double) depth + 4.0 + (double) C.b / 16.0);
            flint_printf("   m %7wd b %4wd np %wd : %10.1f us "
                         "(%.3f ns/unit)\n", C.m, C.b, C.np, t,
                         1000.0 * t / samp[nsamp].units);
            nsamp++;

            sd_fft_mpn_mulmod_2expp1_scratch_clear(&S);
            sd_fft_mpn_mulmod_2expp1_ctx_clear(&C);
            flint_aligned_free(F);
            flint_free(x);
            flint_free(y);
            flint_free(z);
        }
    }

    /* median per prime count */
    for (np = 2; np <= 8; np++)
    {
        double v[256];
        slong n = 0;

        for (i = 0; i < nsamp; i++)
            if (samp[i].np == np)
                v[n++] = samp[i].t / samp[i].units;
        cnt[np] = n;
        if (n == 0)
        {
            med[np] = 0.0;
            continue;
        }
        for (i = 0; i < n; i++)
            for (j = i + 1; j < n; j++)
                if (v[j] < v[i])
                {
                    double s = v[i];
                    v[i] = v[j];
                    v[j] = s;
                }
        med[np] = v[n / 2];
    }

    /* baseline: pooled np = 3..4 median (fall back to any measured) */
    base = 0.0;
    if (cnt[3] || cnt[4])
        base = (med[3] * cnt[3] + med[4] * cnt[4])
               / (double) (cnt[3] + cnt[4]);
    else
        for (np = 2; np <= 8; np++)
            if (cnt[np] && (base == 0.0 || med[np] < base))
                base = med[np];

    flint_printf("\nper-unit medians (ns):");
    for (np = 2; np <= 8; np++)
        if (cnt[np])
            flint_printf("  np%wd %.3f (x%wd)", np,
                         1000.0 * med[np], cnt[np]);
    flint_printf("\n\nsuggested weights (paste into "
                 "_mulmod_2expp1_choose_m, src/fft_small/negacyclic.c;\n"
                 "unmeasured prime counts inherit 14 and are marked):\n\n");
    flint_printf("    static const slong npw[9] = { 0, 0");
    for (np = 2; np <= 8; np++)
    {
        w[np] = cnt[np]
              ? (slong) (16.0 * med[np] / base * 14.0 / 16.0 + 0.5)
              : 14;
        if (w[np] < 8)
            w[np] = 8;
        if (w[np] > 64)
            w[np] = 64;
        flint_printf(", %wd%s", w[np], cnt[np] ? "" : " /*unmeasured*/");
    }
    flint_printf(" };\n\n");

    /* winner check */
    for (i = 0; i < nN; i++)
    {
        slong N = Ns[i], bestt = -1, bestw = -1;
        double tt = 0, tw = 0;

        for (j = 0; j < nsamp; j++)
        {
            double sc;

            if (samp[j].N != N)
                continue;
            sc = (double) w[samp[j].np] * samp[j].units;
            if (bestt < 0 || samp[j].t < tt)
            {
                bestt = j;
                tt = samp[j].t;
            }
            if (bestw < 0 || sc < tw)
            {
                bestw = j;
                tw = sc;
            }
        }
        if (bestt < 0)
            continue;
        flint_printf("N %8wd: measured winner b %3wd np %wd, "
                     "weighted choice b %3wd np %wd  %s",
                     N, samp[bestt].b, samp[bestt].np,
                     samp[bestw].b, samp[bestw].np,
                     bestt == bestw ? "[agree]\n" : "");
        if (bestt != bestw)
            flint_printf("[DISAGREE: %.1f vs %.1f us, %+.1f%%]\n",
                         samp[bestw].t, samp[bestt].t,
                         100.0 * (samp[bestw].t / samp[bestt].t - 1.0));
    }

    flint_cleanup_master();
    return 0;
}
