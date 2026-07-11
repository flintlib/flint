/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "flint.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "profiler.h"

typedef void (*mulmid_fn)(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t,
                          mp_size_t, mp_size_t);

#define IDX_DISPATCH 0
#define NALG 7

static const char * const alg_name[NALG] =
{
    "dispatch", "classical", "via_mul", "via_mullow",
    "via_mulhigh", "via_n_pad", "fft_small"
};

static const mulmid_fn alg_fn[NALG] =
{
    flint_mpn_mulmid,
    flint_mpn_mulmid_classical,
    flint_mpn_mulmid_via_mul,
    flint_mpn_mulmid_via_mullow_n,
    flint_mpn_mulmid_via_mulhigh_n,
    flint_mpn_mulmid_via_n_padded,
    flint_mpn_mulmid_fft_small
};

/* seconds per call, adaptively repeated until the cpu clock is meaningful */
static double
time_call(mulmid_fn f, mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b,
          mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
{
    timeit_t timer;
    slong reps;

    f(z, a, an, b, bn, zlo, zhi);       /* warm up / prime caches */

    TIMEIT_REPEAT(timer, reps)
        f(z, a, an, b, bn, zlo, zhi);
    TIMEIT_END_REPEAT(timer, reps);

    return (double) timer->wall * 0.001 / reps;
}

/* shape categories, chosen to bracket the algorithm crossovers */
enum {
    SH_LOW, SH_NEAR_LOW,
    SH_HIGH, SH_NEAR_HIGH, SH_BAL_MID, SH_NEAR_BAL_MID,
    SH_RANDOM, SH_NEAR_FULL,
    SH_NCAT
};

static const char * const cat_name[SH_NCAT] =
{
    "low", "near-low",
    "high", "near-high", "bal-mid", "near-bal-mid",
    "near-full",
    "random",
};

/* small signed jitter in [-j, j] */
static mp_size_t
jitter(flint_rand_t state, mp_size_t j)
{
    return (mp_size_t) (n_randint(state, 2 * j + 1)) - j;
}

static mp_size_t
clampsz(mp_size_t x, mp_size_t lo, mp_size_t hi)
{
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return x;
}

/* fill (an, bn, zlo, zhi) for a category at a given scale */
static void
gen_shape(flint_rand_t state, int cat, mp_size_t scale,
          mp_size_t * an, mp_size_t * bn, mp_size_t * zlo, mp_size_t * zhi)
{
    mp_size_t a, b, n, M, tot, lo, hi;

    switch (cat)
    {
    case SH_RANDOM:
        a = 1 + n_randint(state, scale);
        b = 1 + n_randint(state, scale);
        tot = a + b;
        lo = n_randint(state, tot);
        hi = lo + 1 + n_randint(state, tot - lo);
        break;

    case SH_NEAR_FULL:                       /* whole product but a few limbs */
        a = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        b = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        tot = a + b;
        lo = n_randint(state, 5);
        hi = tot - n_randint(state, 5);
        break;

    case SH_LOW:                             /* exact low product, mullow shape */
        a = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        b = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        tot = a + b;
        lo = 0;
        hi = clampsz(FLINT_MIN(a, b) + jitter(state, 4), 1, tot);
        break;

    case SH_NEAR_LOW:                        /* low but with a small nonzero zlo */
        a = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        b = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        tot = a + b;
        lo = 1 + n_randint(state, 4);
        hi = clampsz(FLINT_MIN(a, b) + jitter(state, 4), lo + 1, tot);
        break;

    case SH_HIGH:                            /* top slice, mulhigh shape */
        a = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        b = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        tot = a + b;
        M = FLINT_MAX(a, b);
        lo = M + n_randint(state, FLINT_MIN(a, b));
        hi = tot;
        break;

    case SH_NEAR_HIGH:                       /* near the mulhigh boundary */
        a = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        b = scale / 2 + 1 + n_randint(state, scale / 2 + 1);
        tot = a + b;
        M = FLINT_MAX(a, b);
        lo = clampsz(M + jitter(state, 4), 0, tot - 2);
        hi = tot - n_randint(state, 5);
        break;

    case SH_BAL_MID:                         /* balanced middle: mulmid_n shape */
        n = FLINT_MAX((mp_size_t) 2, scale);
        a = 2 * n - 1;
        b = n;
        lo = n - 1;
        hi = 2 * n - 1;
        break;

    case SH_NEAR_BAL_MID:                    /* perturbed balanced middle */
    default:
        n = FLINT_MAX((mp_size_t) 8, scale);
        a = 2 * n - 1 + jitter(state, 5);
        b = n + jitter(state, 5);
        tot = a + b;
        lo = clampsz(n - 1 + jitter(state, 5), 0, tot - 2);
        hi = clampsz(2 * n - 1 + jitter(state, 5), lo + 1, tot);
        break;
    }

    if (hi > a + b) hi = a + b;
    if (lo >= hi)   lo = hi - 1;
    if (lo < 0)     lo = 0;

    *an = a; *bn = b; *zlo = lo; *zhi = hi;
}

int main(int argc, char ** argv)
{
    flint_rand_t state;
    mp_size_t scales[] = { 4, 8, 16, 32, 64, 160, 400, 1000, 2500, 10000, };
    slong nscales = sizeof(scales) / sizeof(scales[0]);
    slong reps_per_cell = (argc > 1) ? atol(argv[1]) : 2;
    mp_size_t maxn = 2 * scales[nscales - 1] + 16;
    mp_ptr a, b, z, zref, d;
    double worst_ratio = 1.0;
    char worst_desc[128];
    int cat;
    slong s, r;

    flint_rand_init(state);

    a    = flint_malloc(sizeof(mp_limb_t) * maxn);
    b    = flint_malloc(sizeof(mp_limb_t) * maxn);
    z    = flint_malloc(sizeof(mp_limb_t) * 2 * maxn);
    zref = flint_malloc(sizeof(mp_limb_t) * 2 * maxn);
    d    = flint_malloc(sizeof(mp_limb_t) * 2 * maxn);
    strcpy(worst_desc, "(none)");

    for (cat = 0; cat < SH_NCAT; cat++)
    {
        flint_printf("%-14s %6s %6s %8s %8s |", "shape", "an", "bn", "zlo", "zn");
        for (int j = 0; j < NALG; j++)
            flint_printf(" %10s", alg_name[j]);
        flint_printf(" | %8s %s\n", "disp/best", "best");

        for (s = 0; s < nscales; s++)
        for (r = 0; r < reps_per_cell; r++)
        {
            mp_size_t an, bn, zlo, zhi, zn;
            double t[NALG];
            double best = 0.0;
            int best_j = 0, j;

            do {
                gen_shape(state, cat, scales[s], &an, &bn, &zlo, &zhi);
            } while (zhi == 0);

            zn = zhi - zlo;

            flint_mpn_rrandom(a, state, an);
            flint_mpn_rrandom(b, state, bn);

            /* exact reference (high limbs) to catch gross errors while tuning */
            flint_mpn_mulmid_via_mul(zref, a, an, b, bn, zlo, zhi);

            for (j = 0; j < NALG; j++)
            {
                t[j] = time_call(alg_fn[j], z, a, an, b, bn, zlo, zhi);

                /* deficit = exact - computed must be non-negative and confined to
                   the low two window limbs (borrow through zero runs may still make
                   individual higher limbs differ, so compare the deficit, not the
                   limbs) */
                {
                    int bad = 0;
                    mp_size_t i;
                    mpn_sub_n(d, zref, z, zn);   /* deficit mod 2^(64*zn) */
                    for (i = 2; i < zn && !bad; i++)
                        if (d[i] != 0)
                            bad = 1;
                    if (zn >= 2 && d[1] > (mp_limb_t) (an + bn))
                        bad = 1;
                    if (bad)
                        flint_printf("MISMATCH %s an=%wd bn=%wd zlo=%wd zhi=%wd\n",
                                     alg_name[j], an, bn, zlo, zhi);
                }

                if (j == 0 || t[j] < best) { best = t[j]; best_j = j; }
            }

            {
                double ratio = (best > 0.0) ? t[IDX_DISPATCH] / best : 1.0;

                flint_printf("%-14s %6wd %6wd %8wd %8wd |",
                             cat_name[cat], an, bn, zlo, zn);
                for (j = 0; j < NALG; j++)
                    flint_printf(" %10.2e", t[j]);
                flint_printf(" | %8.2f %s\n", ratio, alg_name[best_j]);

                if (ratio > worst_ratio && best_j != IDX_DISPATCH)
                {
                    worst_ratio = ratio;
                    flint_sprintf(worst_desc, "%s an=%wd bn=%wd zlo=%wd zn=%wd (best %s)",
                                  cat_name[cat], an, bn, zlo, zn, alg_name[best_j]);
                }
            }
        }
    }

    flint_printf("\nworst dispatch slowdown: %.2fx at %s\n", worst_ratio, worst_desc);

    flint_free(a); flint_free(b); flint_free(z); flint_free(zref); flint_free(d);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
