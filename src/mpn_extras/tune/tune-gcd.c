/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Tuning program for the gcd and half-gcd cutoffs in flint-mparam.h:

        FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF   classical -> Strassen-like
                                            2 x 2 matrix products
        FLINT_MPN_HGCD_CUTOFF               Lehmer -> recursive half-gcd
        FLINT_MPN_HGCD_APPR_CUTOFF          the same, approximate half-gcd
        FLINT_MPN_GCD_DC_CUTOFF             Lehmer -> half-gcd based gcd
        FLINT_MPN_GCDEXT_DC_CUTOFF          the same, extended gcd (both by
                                            default derived from
                                            FLINT_MPN_HGCD_CUTOFF)
        FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF   matrix operations of the
        FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF       half-gcd in fft_small
        FLINT_MPN_HGCD_APPLY_FFT_CUTOFF        transform space (one value,
                                            chosen for all three)
        FLINT_MPN_HGCD_REDUCE_CUTOFF        half-gcd of the top half ->
                                            approximate half-gcd in the
                                            first recursive step
        FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF  threads in the transform-space
                                            operations (only when run with
                                            more than one thread)

    The gcd sources are compiled into this program with the cutoffs turned
    into variables, so that all crossovers are measured in one run without
    rebuilding FLINT. Build and run it from the FLINT build directory with

        make build/mpn_extras/tune/tune-gcd
        build/mpn_extras/tune/tune-gcd [threads]

    The cutoffs other than the parallel one are measured with one thread;
    with an argument threads > 1, the parallel cutoff is then measured with
    that many threads. Run on an idle machine; a run takes about half a
    minute and a half, plus about a minute for the parallel cutoff.
    The crossovers are measured in the order listed, each with the values
    found before it. The suggested values are printed at the end. The
    Strassen and approximate half-gcd cutoffs are one-level crossovers,
    whose tables show the time ratios (new method / old method, below 1
    means the new method is faster). The other cutoffs, whose effect is
    only seen on the whole algorithm or whose one-level crossover is flat,
    are swept over candidate values on the gcd and the extended gcd; their
    tables show the time ratios against not using them and their geometric
    mean over the sizes shown, the value the defaults give being marked
    with *, and they are only suggested when the best value is more than
    1% faster than the default. The half-gcd cutoff is swept with the gcd
    cutoffs following it, as their defaults do.
*/

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "flint.h"
#include "thread_support.h"

#if FLINT_BITS != 64
# error "the tuned code is only used on 64-bit machines"
#endif

/* "off": larger than any size used, small enough that the scratch size
   formulas, which add some of the cutoffs, do not overflow */
#define TUNE_OFF ((slong) 1 << 40)

/* cutoffs as variables, initialised from the defaults below */
static slong tune_strassen, tune_hgcd, tune_hgcd_appr, tune_hgcd_reduce,
    tune_gcd_dc, tune_gcdext_dc, tune_fft, tune_parallel;

#include "mpn_extras.h"

/* whether flint-mparam.h sets the gcd cutoffs, which otherwise follow
   from the half-gcd cutoff */
#ifdef FLINT_MPN_GCD_DC_CUTOFF
# define TUNE_GCD_DC_FIXED 1
#else
# define TUNE_GCD_DC_FIXED 0
#endif
#ifdef FLINT_MPN_GCDEXT_DC_CUTOFF
# define TUNE_GCDEXT_DC_FIXED 1
#else
# define TUNE_GCDEXT_DC_FIXED 0
#endif

/* the defaults (flint-mparam.h and impl.h; the transform-space ones as
   in hgcd_fft.c) */
#include "../impl.h"
static const slong default_strassen = FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF;
static const slong default_hgcd = FLINT_MPN_HGCD_CUTOFF;
static const slong default_hgcd_appr = FLINT_MPN_HGCD_APPR_CUTOFF;
static const slong default_hgcd_reduce = FLINT_MPN_HGCD_REDUCE_CUTOFF;
static const slong default_gcd_dc = FLINT_MPN_GCD_DC_CUTOFF;
static const slong default_gcdext_dc = FLINT_MPN_GCDEXT_DC_CUTOFF;
#ifdef FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF
static const slong default_fft = FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF;
#else
static const slong default_fft = FLINT_MAX(150, 3 * FLINT_FFT_SMALL_MUL_THRESHOLD / 10);
#endif
#ifdef FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF
static const slong default_parallel = FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF;
#else
static const slong default_parallel = 500;
#endif

#undef FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF
#undef FLINT_MPN_HGCD_CUTOFF
#undef FLINT_MPN_HGCD_APPR_CUTOFF
#undef FLINT_MPN_HGCD_REDUCE_CUTOFF
#undef FLINT_MPN_GCD_DC_CUTOFF
#undef FLINT_MPN_GCDEXT_DC_CUTOFF
#undef FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF
#undef FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF
#undef FLINT_MPN_HGCD_APPLY_FFT_CUTOFF
#undef FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF
#define FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF tune_strassen
#define FLINT_MPN_HGCD_CUTOFF tune_hgcd
#define FLINT_MPN_HGCD_APPR_CUTOFF tune_hgcd_appr
#define FLINT_MPN_HGCD_REDUCE_CUTOFF tune_hgcd_reduce
#define FLINT_MPN_GCD_DC_CUTOFF tune_gcd_dc
#define FLINT_MPN_GCDEXT_DC_CUTOFF tune_gcdext_dc
#define FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF tune_fft
#define FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF tune_fft
#define FLINT_MPN_HGCD_APPLY_FFT_CUTOFF tune_fft
#define FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF tune_parallel

/* private copies of the exported functions */
#define _flint_mpn_hgcd_matrix_init tune_mpn_hgcd_matrix_init
#define _flint_mpn_hgcd_matrix_update_q tune_mpn_hgcd_matrix_update_q
#define _flint_mpn_hgcd_matrix_mul_1 tune_mpn_hgcd_matrix_mul_1
#define _flint_mpn_hgcd_matrix22_mul_strassen_itch tune_mpn_hgcd_matrix22_mul_strassen_itch
#define _flint_mpn_hgcd_matrix22_mul_strassen tune_mpn_hgcd_matrix22_mul_strassen
#define _flint_mpn_hgcd_matrix_mul tune_mpn_hgcd_matrix_mul
#define _flint_mpn_hgcd_matrix_adjust tune_mpn_hgcd_matrix_adjust
#define _flint_mpn_gcd_subdiv_step tune_mpn_gcd_subdiv_step
#define _flint_mpn_hgcd_step tune_mpn_hgcd_step
#define _flint_mpn_hgcd_itch tune_mpn_hgcd_itch
#define _flint_mpn_hgcd tune_mpn_hgcd
#define _flint_mpn_hgcd_appr_itch tune_mpn_hgcd_appr_itch
#define _flint_mpn_hgcd_appr tune_mpn_hgcd_appr
#define _flint_mpn_hgcd_reduce_itch tune_mpn_hgcd_reduce_itch
#define _flint_mpn_hgcd_reduce tune_mpn_hgcd_reduce
#define _flint_mpn_hgcd_matrix_mul_fft tune_mpn_hgcd_matrix_mul_fft
#define _flint_mpn_hgcd_matrix_adjust_fft tune_mpn_hgcd_matrix_adjust_fft
#define _flint_mpn_hgcd_matrix_adjust_mul_fft tune_mpn_hgcd_matrix_adjust_mul_fft
#define _flint_mpn_hgcd_mul_matrix_vector_fft tune_mpn_hgcd_mul_matrix_vector_fft
#define _flint_mpn_hgcd_matrix_apply_fft tune_mpn_hgcd_matrix_apply_fft
#define _flint_mpn_gcd_inplace tune_mpn_gcd_inplace
#define _flint_mpn_gcdext_hook tune_mpn_gcdext_hook
#define _flint_mpn_gcdext_lehmer_n tune_mpn_gcdext_lehmer_n
#define _flint_mpn_gcdext_inplace tune_mpn_gcdext_inplace

/* prototypes for the renamed functions */
void tune_mpn_hgcd_matrix_init(flint_mpn_hgcd_matrix_struct *, mp_size_t, mp_ptr);
void tune_mpn_hgcd_matrix_update_q(flint_mpn_hgcd_matrix_struct *, mp_srcptr, mp_size_t, unsigned int, mp_ptr);
void tune_mpn_hgcd_matrix_mul_1(flint_mpn_hgcd_matrix_struct *, const flint_mpn_hgcd_matrix1_struct *, mp_ptr);
mp_size_t tune_mpn_hgcd_matrix22_mul_strassen_itch(mp_size_t, mp_size_t);
void tune_mpn_hgcd_matrix22_mul_strassen(mp_ptr, mp_ptr, mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_srcptr, mp_srcptr, mp_srcptr, mp_size_t, mp_ptr);
void tune_mpn_hgcd_matrix_mul(flint_mpn_hgcd_matrix_struct *, const flint_mpn_hgcd_matrix_struct *, mp_ptr);
mp_size_t tune_mpn_hgcd_matrix_adjust(const flint_mpn_hgcd_matrix_struct *, mp_size_t, mp_ptr, mp_ptr, mp_size_t, mp_ptr);
mp_size_t tune_mpn_gcd_subdiv_step(mp_ptr, mp_ptr, mp_size_t, mp_size_t, _flint_mpn_gcd_subdiv_step_hook *, void *, mp_ptr);
mp_size_t tune_mpn_hgcd_step(mp_size_t, mp_ptr, mp_ptr, mp_size_t, flint_mpn_hgcd_matrix_struct *, mp_ptr);
mp_size_t tune_mpn_hgcd_itch(mp_size_t);
mp_size_t tune_mpn_hgcd(mp_ptr, mp_ptr, mp_size_t, flint_mpn_hgcd_matrix_struct *, mp_ptr);
mp_size_t tune_mpn_hgcd_appr_itch(mp_size_t);
int tune_mpn_hgcd_appr(mp_ptr, mp_ptr, mp_size_t, flint_mpn_hgcd_matrix_struct *, mp_ptr);
mp_size_t tune_mpn_hgcd_reduce_itch(mp_size_t, mp_size_t);
mp_size_t tune_mpn_hgcd_reduce(flint_mpn_hgcd_matrix_struct *, mp_ptr, mp_ptr, mp_size_t, mp_size_t, mp_ptr);
int tune_mpn_hgcd_matrix_mul_fft(flint_mpn_hgcd_matrix_struct *, const flint_mpn_hgcd_matrix_struct *);
mp_size_t tune_mpn_hgcd_matrix_adjust_fft(const flint_mpn_hgcd_matrix_struct *, mp_size_t, mp_ptr, mp_ptr, mp_size_t);
mp_size_t tune_mpn_hgcd_matrix_adjust_mul_fft(flint_mpn_hgcd_matrix_struct *, const flint_mpn_hgcd_matrix_struct *, mp_size_t, mp_ptr, mp_ptr, mp_size_t);
mp_size_t tune_mpn_hgcd_mul_matrix_vector_fft(const flint_mpn_hgcd_matrix_struct *, mp_ptr, mp_srcptr, mp_ptr, mp_size_t);
int tune_mpn_hgcd_matrix_apply_fft(const flint_mpn_hgcd_matrix_struct *, mp_size_t [2][2], mp_ptr, mp_ptr, mp_size_t, mp_size_t);
mp_size_t tune_mpn_gcd_inplace(mp_ptr, mp_ptr, mp_size_t, mp_ptr, mp_size_t);
void tune_mpn_gcdext_hook(void *, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, int);
mp_size_t tune_mpn_gcdext_lehmer_n(mp_ptr, mp_ptr, mp_size_t *, mp_ptr, mp_ptr, mp_size_t, mp_ptr);
mp_size_t tune_mpn_gcdext_inplace(mp_ptr, mp_ptr, mp_size_t *, mp_ptr, mp_size_t, mp_ptr, mp_size_t);

/* the sources, which see the renamed cutoffs and functions (impl.h and
   mpn_extras.h are not read again) */
#include "../hgcd.c"
#include "../hgcd_fft.c"
#include "../gcd_inplace.c"
#include "../gcdext_inplace.c"

/* timing */
static double
now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

/* time per call of expr1 and expr2, measured alternately (the minimum over
   7 rounds each), to reduce the effect of frequency changes and noise */
#define TIME2(expr1, res1, expr2, res2) \
    do { \
        double __b1 = 1e30, __b2 = 1e30, __t, __t0; \
        long __reps, __i; \
        int __r; \
        for (__reps = 1; ; __reps *= 2) \
        { \
            __t0 = now(); \
            for (__i = 0; __i < __reps; __i++) { expr1; } \
            __t = now() - __t0; \
            if (__t > 0.003) break; \
        } \
        for (__r = 0; __r < 7; __r++) \
        { \
            __t0 = now(); \
            for (__i = 0; __i < __reps; __i++) { expr1; } \
            __t = (now() - __t0) / __reps; \
            if (__t < __b1) __b1 = __t; \
            __t0 = now(); \
            for (__i = 0; __i < __reps; __i++) { expr2; } \
            __t = (now() - __t0) / __reps; \
            if (__t < __b2) __b2 = __t; \
        } \
        (res1) = __b1; \
        (res2) = __b2; \
    } while (0)

static flint_rand_t state;

/* while the half-gcd cutoff is swept, the gcd cutoffs follow it as their
   defaults do (unless flint-mparam.h sets them) */
static int tie_dc = 0;

static void
sync_dc(void)
{
    if (tie_dc)
    {
        if (!TUNE_GCD_DC_FIXED)
            tune_gcd_dc = 3 * tune_hgcd;
        if (!TUNE_GCDEXT_DC_FIXED)
            tune_gcdext_dc = 3 * tune_hgcd / 2;
    }
}

/* operands and scratch, for sizes up to MAXN limbs */
#define MAXN 200000
static mp_ptr A, B, A2, B2, G, S, T, MP;

static void
random_operands(mp_size_t n)
{
    flint_mpn_urandomb(A, state, n * FLINT_BITS);
    flint_mpn_urandomb(B, state, n * FLINT_BITS);
    A[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    B[n - 1] |= UWORD(1) << (FLINT_BITS - 2);
    A[0] |= 1;
}

static void
run_gcd(mp_size_t n)
{
    flint_mpn_copyi(A2, A, n);
    flint_mpn_copyi(B2, B, n);
    sync_dc();
    tune_mpn_gcd_inplace(G, A2, n, B2, n);
}

static void
run_gcdext(mp_size_t n)
{
    mp_size_t sn;
    flint_mpn_copyi(A2, A, n);
    flint_mpn_copyi(B2, B, n);
    sync_dc();
    tune_mpn_gcdext_inplace(G, S, &sn, A2, n, B2, n);
}

static void
run_hgcd(mp_size_t n)
{
    flint_mpn_hgcd_matrix_struct M;
    flint_mpn_copyi(A2, A, n);
    flint_mpn_copyi(B2, B, n);
    tune_mpn_hgcd_matrix_init(&M, n, MP);
    tune_mpn_hgcd(A2, B2, n, &M, T);
}

static void
run_hgcd_appr(mp_size_t n)
{
    flint_mpn_hgcd_matrix_struct M;
    flint_mpn_copyi(A2, A, n);
    flint_mpn_copyi(B2, B, n);
    tune_mpn_hgcd_matrix_init(&M, n, MP);
    tune_mpn_hgcd_appr(A2, B2, n, &M, T);
}

/*
    Crossover from a table of ratios (new / old): the first size from which
    the new method is faster and stays within 2% of the old one at all larger
    sizes, ignoring isolated outliers (a single size whose neighbours are
    both faster). Returns -1 if there is none.
*/
static slong
crossover(const slong * sizes, const double * ratio, int num)
{
    int i, j, ok;

    for (i = 0; i < num; i++)
    {
        if (ratio[i] >= 1.0)
            continue;
        ok = 1;
        for (j = i + 1; j < num; j++)
        {
            if (ratio[j] > 1.02
                && !(j + 1 < num && ratio[j - 1] < 1.0 && ratio[j + 1] < 1.0))
                ok = 0;
        }
        if (ok)
            return sizes[i];
    }
    return -1;
}

static void
print_row(slong size, double t_old, double t_new)
{
    flint_printf("  %6wd   %10.3e %10.3e   %6.3f\n", size, t_old, t_new, t_new / t_old);
}

#define MAXS 64

/* One-level crossover of a cutoff: the time of run(n) with the cutoff
   off against the cutoff n (so that only the top level uses the new
   method), for the given sizes; returns the crossover, or fallback. */
static slong
tune_one_level(const char * name, const char * what, slong * var,
    void (*run)(mp_size_t), slong n0, slong n1, double step, slong fallback)
{
    slong sizes[MAXS], n, result;
    double ratio[MAXS], t1, t2;
    int num = 0;

    flint_printf("%s: %s\n", name, what);
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "old", "new", "ratio");

    for (n = n0; n <= n1 && num < MAXS; n = FLINT_MAX(n + 1, (slong) (n * step)))
    {
        random_operands(n);
        TIME2((*var = TUNE_OFF, run(n)), t1, (*var = n, run(n)), t2);
        print_row(n, t1, t2);
        sizes[num] = n;
        ratio[num++] = t2 / t1;
    }

    result = crossover(sizes, ratio, num);
    result = (result == -1) ? fallback : result;
    *var = result;
    flint_printf("  -> %wd\n\n", result);
    return result;
}

/* Sweep of a cutoff over the candidate values (the last one being
   TUNE_OFF) and the reference value ref (the value the defaults give),
   timing the gcd (which & 1) and the extended gcd (which & 2) at the
   given sizes, all candidates measured in turn, the minimum over 5
   rounds. Returns the candidate with the smallest geometric mean of the
   times relative to TUNE_OFF, and sets *ref_excess to the ratio of the
   mean of ref to that of the best candidate. */
static slong
tune_sweep(const char * name, slong * var, const slong * cand0, int ncand0,
    slong ref, const slong * sizes, int nsizes, int which, double * ref_excess)
{
    double t[20][MAXS][2], logsum[20], noise;
    slong cand[20];
    int c, k, f, r, best, iref, ncand, nnoise;
    long reps, i;
    double t0, tt;

    /* the candidates with ref inserted in order, TUNE_OFF last */
    ncand = 0;
    for (c = 0; c < ncand0 - 1; c++)
    {
        if (ref < cand0[c] && (ncand == 0 || cand[ncand - 1] < ref))
            cand[ncand++] = ref;
        if (cand0[c] != ref)
            cand[ncand++] = cand0[c];
        else if (ncand == 0 || cand[ncand - 1] != ref)
            cand[ncand++] = ref;
    }
    if (ref < TUNE_OFF && (ncand == 0 || cand[ncand - 1] < ref))
        cand[ncand++] = ref;
    cand[ncand++] = TUNE_OFF;

    flint_printf("%s: geometric mean of the time ratios against not using it, %s\n",
        name, (which == 3) ? "gcd and gcdext" : (which == 1) ? "gcd" : "gcdext");

    for (k = 0; k < nsizes; k++)
    {
        mp_size_t n = sizes[k];
        random_operands(n);

        for (f = 0; f < 2; f++)
        {
            void (*run)(mp_size_t) = (f == 0) ? run_gcd : run_gcdext;

            for (c = 0; c < ncand; c++)
                t[c][k][f] = 1.0;

            if (!(which & (1 << f)))
                continue;

            for (c = 0; c < ncand; c++)
                t[c][k][f] = 1e30;

            *var = cand[ncand - 1];
            for (reps = 1; ; reps *= 2)
            {
                t0 = now();
                for (i = 0; i < reps; i++)
                    run(n);
                if (now() - t0 > 0.003)
                    break;
            }

            for (r = 0; r < 5; r++)
            {
                for (c = 0; c < ncand; c++)
                {
                    *var = cand[c];
                    t0 = now();
                    for (i = 0; i < reps; i++)
                        run(n);
                    tt = (now() - t0) / reps;
                    if (tt < t[c][k][f])
                        t[c][k][f] = tt;
                }
            }
        }
    }

    flint_printf("  %8s", "cutoff");
    for (k = 0; k < nsizes; k++)
        flint_printf(" %6wd", sizes[k]);
    flint_printf("   %6s\n", "mean");

    iref = ncand - 1;
    noise = 0;
    nnoise = 0;
    for (c = 0; c < ncand; c++)
    {
        logsum[c] = 0;
        if (cand[c] == ref)
            iref = c;
        if (cand[c] >= TUNE_OFF)
            flint_printf("  %8s", "off");
        else
            flint_printf("  %7wd%s", cand[c], (cand[c] == ref) ? "*" : " ");
        for (k = 0; k < nsizes; k++)
        {
            double q = (t[c][k][0] / t[ncand - 1][k][0]) * (t[c][k][1] / t[ncand - 1][k][1]);
            if (which == 3)
                q = sqrt(q);
            flint_printf(" %6.3f", q);
            logsum[c] += log(q);
            /* a cutoff above the size gives the same algorithm as off, so
               these entries measure the noise */
            if (cand[c] > sizes[k] && cand[c] < TUNE_OFF)
            {
                noise += log(q) * log(q);
                nnoise++;
            }
        }
        logsum[c] /= nsizes;
        flint_printf("   %6.3f\n", exp(logsum[c]));
    }

    best = ncand - 1;
    for (c = 0; c < ncand; c++)
        if (logsum[c] < logsum[best])
            best = c;

    *ref_excess = exp(logsum[iref] - logsum[best]);
    /* the following measurements use the value that will be suggested */
    *var = (*ref_excess > 1.01) ? cand[best] : ref;
    if (cand[best] >= TUNE_OFF)
        flint_printf("  -> off");
    else
        flint_printf("  -> %wd", cand[best]);
    flint_printf(" (the default * is %.1f%% slower)\n", 100 * (*ref_excess - 1));
    if (nnoise > 0)
    {
        noise = 100 * (exp(sqrt(noise / nnoise)) - 1);
        flint_printf("  noise: %.1f%% (the entries for cutoffs above the size)%s\n", noise,
            (noise > 2.0) ? "; too noisy, run again on an idle machine" : "");
    }
    flint_printf("\n");
    return cand[best];
}

int main(int argc, char ** argv)
{
    slong sizes[MAXS], n, result, threads, ref_gcd, ref_gcdext, best_hgcd, best_reduce;
    double ratio[MAXS], t1, t2, ex_hgcd, ex_gcd, ex_gcdext, ex_reduce, ex_fft = 1, ex_par = 1;
    int num;

    /* the serial cutoffs are measured with one thread, the parallel cutoff
       with the given number */
    threads = (argc > 1) ? atol(argv[1]) : 1;
    flint_set_num_threads(1);
    flint_rand_init(state);

    A = flint_malloc((MAXN + 1) * sizeof(mp_limb_t));
    B = flint_malloc((MAXN + 1) * sizeof(mp_limb_t));
    A2 = flint_malloc((MAXN + 1) * sizeof(mp_limb_t));
    B2 = flint_malloc((MAXN + 1) * sizeof(mp_limb_t));
    G = flint_malloc((MAXN + 1) * sizeof(mp_limb_t));
    S = flint_malloc((MAXN + 1) * sizeof(mp_limb_t));
    /* scratch for the half-gcd runs */
    T = flint_malloc((30 * 20000 + 10000) * sizeof(mp_limb_t));
    MP = flint_malloc((4 * 20000 + 100) * sizeof(mp_limb_t));

    tune_strassen = default_strassen;
    tune_hgcd = default_hgcd;
    tune_hgcd_appr = default_hgcd_appr;
    tune_hgcd_reduce = default_hgcd_reduce;
    tune_gcd_dc = default_gcd_dc;
    tune_gcdext_dc = default_gcdext_dc;
    tune_fft = default_fft;
    tune_parallel = TUNE_OFF;

    /* 1. the Strassen-like 2 x 2 product */
    flint_printf("FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF: n x n matrix product, classical vs Strassen-like\n");
    flint_printf("  %6s   %10s %10s   %6s\n", "n", "classical", "strassen", "ratio");
    num = 0;
    for (n = 6; n <= 80; n += (n < 30) ? 2 : 6)
    {
        mp_ptr r[4], m[4];
        int i;
        for (i = 0; i < 4; i++)
        {
            r[i] = A2 + i * (2 * n + 2);
            m[i] = B2 + i * n;
            flint_mpn_urandomb(r[i], state, n * FLINT_BITS);
            flint_mpn_urandomb(m[i], state, n * FLINT_BITS);
        }
        TIME2((tune_strassen = TUNE_OFF, tune_mpn_hgcd_matrix22_mul_strassen(r[0], r[1], r[2], r[3], n, m[0], m[1], m[2], m[3], n, T)), t1,
              (tune_strassen = n, tune_mpn_hgcd_matrix22_mul_strassen(r[0], r[1], r[2], r[3], n, m[0], m[1], m[2], m[3], n, T)), t2);
        print_row(n, t1, t2);
        sizes[num] = n; ratio[num++] = t2 / t1;
    }
    result = crossover(sizes, ratio, num);
    tune_strassen = (result == -1) ? 80 : result;
    flint_printf("  -> %wd\n\n", tune_strassen);

    /* 2.-5. Lehmer vs the half-gcd based algorithms. The half-gcd cutoff
       is swept on the gcd and the extended gcd with the gcd cutoffs
       following it, as their defaults do: its one-level crossover is flat
       and moves between runs (77 to 101 on Zen 3). The gcd cutoffs also
       decide where the half-gcd based algorithm hands over to Lehmer's
       (in the extended gcd, with long cofactors), so they are swept on the
       whole algorithm too, to check the values that follow. */
    {
        static const slong cand[] = { 60, 70, 80, 90, 100, 115, 130, 150, 175, 200, 250, TUNE_OFF };
        static const slong hsizes[] = { 150, 200, 260, 340, 450, 600, 800, 1100, 1500, 2000, 3000 };
        tie_dc = 1;
        best_hgcd = tune_sweep("FLINT_MPN_HGCD_CUTOFF (the gcd cutoffs following it)", &tune_hgcd, cand, 12, default_hgcd, hsizes, 11, 3, &ex_hgcd);
        sync_dc();
        tie_dc = 0;
    }
    tune_one_level("FLINT_MPN_HGCD_APPR_CUTOFF", "approximate half-gcd of n limbs, Lehmer vs recursive",
        &tune_hgcd_appr, run_hgcd_appr, 20, 1500, 1.15, 1500);
    ref_gcd = TUNE_GCD_DC_FIXED ? default_gcd_dc : 3 * tune_hgcd;
    ref_gcdext = TUNE_GCDEXT_DC_FIXED ? default_gcdext_dc : 3 * tune_hgcd / 2;
    {
        static const slong cand[] = { 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, TUNE_OFF };
        static const slong dsizes[] = { 200, 260, 340, 450, 600, 800, 1100, 1500, 2000, 3000 };
        tune_sweep("FLINT_MPN_GCD_DC_CUTOFF", &tune_gcd_dc, cand, 12, ref_gcd, dsizes, 10, 1, &ex_gcd);
    }
    {
        static const slong cand[] = { 120, 150, 180, 210, 250, 300, 350, 400, 500, 600, 800, TUNE_OFF };
        static const slong dsizes[] = { 150, 200, 260, 340, 450, 600, 800, 1100, 1500, 2000, 3000 };
        tune_sweep("FLINT_MPN_GCDEXT_DC_CUTOFF", &tune_gcdext_dc, cand, 12, ref_gcdext, dsizes, 11, 2, &ex_gcdext);
    }

    /* 6. the transform-space matrix operations */
#if FLINT_HAVE_FFT_SMALL
    {
        static const slong cand[] = { 40, 50, 64, 80, 100, 125, 150, 200, 250, 300, 400, TUNE_OFF };
        static const slong fsizes[] = { 150, 200, 260, 340, 450, 600, 800, 1100, 1500, 2000, 3000 };
        tune_sweep("FLINT_MPN_HGCD_*_FFT_CUTOFF", &tune_fft, cand, 12, default_fft, fsizes, 11, 3, &ex_fft);
    }
#endif

    /* 7. the reduction in the first recursive step of the half-gcd; its
       one-level table is flat (within 2% from 600 to 10000 limbs), so it is
       swept on the whole algorithm */
    {
        static const slong cand[] = { 600, 800, 1000, 1300, 1700, 2200, 3000, 4000, 6000, TUNE_OFF };
        static const slong rsizes[] = { 1000, 1500, 2200, 3300, 5000, 7500, 11000, 16000 };
        best_reduce = tune_sweep("FLINT_MPN_HGCD_REDUCE_CUTOFF", &tune_hgcd_reduce, cand, 10, default_hgcd_reduce, rsizes, 8, 3, &ex_reduce);
    }

    /* 8. threads in the transform-space operations */
#if FLINT_HAVE_FFT_SMALL
    if (threads > 1)
    {
        static const slong cand[] = { 250, 500, 1000, 2000, 4000, 8000, TUNE_OFF };
        static const slong psizes[] = { 10000, 30000, 100000 };
        flint_set_num_threads(threads);
        tune_sweep("FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF", &tune_parallel, cand, 7, default_parallel, psizes, 3, 3, &ex_par);
    }
    else
        tune_parallel = default_parallel;
#endif

    /* The swept cutoffs are only suggested when the best value is more
       than 1% faster than the default; the one-level crossovers always
       (their tables show whether they matter). */
    flint_printf("Suggested values for flint-mparam.h (defaults in parentheses):\n\n");
    flint_printf("#define FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF %wd   /* (%wd) */\n", tune_strassen, default_strassen);
    if (best_hgcd >= TUNE_OFF)
        best_hgcd = WORD_MAX;
    if (ex_hgcd > 1.01)
        flint_printf("#define FLINT_MPN_HGCD_CUTOFF %wd   /* (%wd) */\n", best_hgcd, default_hgcd);
    else
        flint_printf("/* FLINT_MPN_HGCD_CUTOFF: the default %wd is within 1%% of the best */\n", default_hgcd);
    flint_printf("#define FLINT_MPN_HGCD_APPR_CUTOFF %wd   /* (%wd) */\n", tune_hgcd_appr, default_hgcd_appr);
    if (best_reduce >= TUNE_OFF)
        best_reduce = WORD_MAX;
    if (ex_reduce > 1.01)
        flint_printf("#define FLINT_MPN_HGCD_REDUCE_CUTOFF %wd   /* (%wd) */\n", best_reduce, default_hgcd_reduce);
    else
        flint_printf("/* FLINT_MPN_HGCD_REDUCE_CUTOFF: the default %wd is within 1%% of the best */\n", default_hgcd_reduce);
    if (tune_gcd_dc >= TUNE_OFF)
        tune_gcd_dc = WORD_MAX;
    if (tune_gcdext_dc >= TUNE_OFF)
        tune_gcdext_dc = WORD_MAX;
    if (ex_gcd > 1.01)
        flint_printf("#define FLINT_MPN_GCD_DC_CUTOFF %wd   /* (%wd) */\n", tune_gcd_dc, ref_gcd);
    else
        flint_printf("/* FLINT_MPN_GCD_DC_CUTOFF: the default %wd%s is within 1%% of the best */\n",
            ref_gcd, TUNE_GCD_DC_FIXED ? "" : " (3 FLINT_MPN_HGCD_CUTOFF)");
    if (ex_gcdext > 1.01)
        flint_printf("#define FLINT_MPN_GCDEXT_DC_CUTOFF %wd   /* (%wd) */\n", tune_gcdext_dc, ref_gcdext);
    else
        flint_printf("/* FLINT_MPN_GCDEXT_DC_CUTOFF: the default %wd%s is within 1%% of the best */\n",
            ref_gcdext, TUNE_GCDEXT_DC_FIXED ? "" : " (3 FLINT_MPN_HGCD_CUTOFF / 2)");
#if FLINT_HAVE_FFT_SMALL
    if (tune_fft >= TUNE_OFF)
        tune_fft = WORD_MAX;
    if (ex_fft > 1.01)
    {
        flint_printf("#define FLINT_MPN_HGCD_MATRIX_MUL_FFT_CUTOFF %wd   /* (%wd) */\n", tune_fft, default_fft);
        flint_printf("#define FLINT_MPN_HGCD_ADJUST_FFT_CUTOFF %wd   /* (%wd) */\n", tune_fft, default_fft);
        flint_printf("#define FLINT_MPN_HGCD_APPLY_FFT_CUTOFF %wd   /* (%wd) */\n", tune_fft, default_fft);
    }
    else
        flint_printf("/* FLINT_MPN_HGCD_*_FFT_CUTOFF: the default %wd is within 1%% of the best */\n", default_fft);
    if (threads > 1)
    {
        if (tune_parallel >= TUNE_OFF)
            tune_parallel = WORD_MAX;
        if (ex_par > 1.01)
            flint_printf("#define FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF %wd   /* (%wd) */\n", tune_parallel, default_parallel);
        else
            flint_printf("/* FLINT_MPN_HGCD_FFT_PARALLEL_CUTOFF: the default %wd is within 1%% of the best */\n", default_parallel);
    }
#endif

    flint_free(A);
    flint_free(B);
    flint_free(A2);
    flint_free(B2);
    flint_free(G);
    flint_free(S);
    flint_free(T);
    flint_free(MP);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
