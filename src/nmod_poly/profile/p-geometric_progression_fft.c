/*
   Copyright (C) 2026 Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
   This can be used to get a good overall view of the stored-FFT impact:
P="taskset -c 2 build/nmod_poly/profile/p-geometric_progression_fft"
$P                              > table.txt       # 6 moduli, full output
$P -s 4096 20 fft ntt 60        > sweep4096.txt   # no-slowdown check
$P -s 3000 fft 60               > sweep3000.txt   # G length not a power of two

   To eliminate disturbing cache / page faults effects from benchmarks
   (especially in the progression init timings), you might use e.g.
GLIBC_TUNABLES=glibc.malloc.trim_threshold=1000000000:glibc.malloc.mmap_threshold=1000000000 \
  taskset -c 2 build/nmod_poly/profile/p-geometric_progression_fft > table.txt
*/

/*
   Geometric evaluation, interpolation and extrapolation with a precomputation
   G built without and with (bit 3 of function) stored transforms.

   usage: p-geometric_progression_fft [--concise] [-s L] [mod ...]

   where each mod is a bit length (a random prime of that size is used),
   "fft" (a prime of the default fft_small context, for which fft_small
   needs a single prime), or "ntt" (998244353 = 119*2^23 + 1, which
   fft_small uses directly, without chinese remaindering).
   Operations for a G of length n: evaluation of a polynomial of length n at
   n points, interpolation from n points, extrapolation from n/2 values to
   the next n/2 points (the largest calls that G allows).

   Without -s, for each length n:
   - times of one call without / with transforms (us) and their ratio, for
     each operation (G with functions 7 and 15);
   - for each precomputation, the times of its init without / with
     transforms and b/e, the number of rounds of one call of each of its
     operations after which the extra init cost is recovered ("-" if no
     transform is stored, "inf" if they do not gain). Precomputations, all
     of length n: for evaluation (function 1), interpolation (2),
     extrapolation (4), evaluation and interpolation (3), all (7).
   With --concise, a single table: the three ratios, and for each
   precomputation the ratio of the init times with / without transforms
   and b/e.

   With -s L: G of length L, and calls for all lengths n <= L on a grid
   (step about n/32), with the shapes above for a G of length n. Reports,
   for each operation, the smallest and largest ratio, how many lengths have
   ratio < 0.97 (these are timed again, longer, before being counted and
   printed), and the largest loss in absolute time.
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

/* minimum over 3 runs of at least ms milliseconds of the time per call, in us */
#define TIME_MIN(t, ms, stmt)                                             \
    do {                                                                  \
        timeit_t __T;                                                     \
        slong __reps = 1;                                                 \
        (t) = 1e300;                                                      \
        for (int __r = 0; __r < 3; )                                      \
        {                                                                 \
            timeit_start_us(__T);                                         \
            for (slong __k = 0; __k < __reps; __k++)                      \
            {                                                             \
                stmt;                                                     \
            }                                                             \
            timeit_stop_us(__T);                                          \
            if (__T->wall < 1000 * (ms))                                  \
                __reps *= 2;                                              \
            else                                                          \
            {                                                             \
                (t) = FLINT_MIN((t), (double) __T->wall / __reps);        \
                __r++;                                                    \
            }                                                             \
        }                                                                 \
    } while (0)

static const slong lens[] = {16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512,
    700, 1024, 1400, 2048, 3000, 4096, 6000, 8192, 12000, 16384};

static const char * names[] = {"eval", "interp", "extrap"};

/* precomputations: name and function (bit i: operation i) */
#define NPRE 5
static const struct
{
    const char * name;
    ulong f;
}
pre[NPRE] = {
    {"eval", 1},
    {"interp", 2},
    {"extrap", 4},
    {"ev+int", 3},
    {"all", 7},
};

/* operation op (0: evaluation, 1: interpolation, 2: extrapolation) for
   length n, as described above */
static void
run_op(int op, nn_ptr b, nn_srcptr a, slong n, const nmod_geometric_progression_struct * G)
{
    slong h = n / 2;

    if (op == 0)
        _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(b, a, n, G, n, G->mod);
    else if (op == 1)
        _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(b, a, G, n, G->mod);
    else
        nmod_poly_extrapolate_geometric_precomp(b, h, a, h, h, G);
}

/* times in us of op for length n with G0 and G1: minimum over 3 alternating
   runs of the same number of calls, at least ms milliseconds with G0 */
static void
time_pair(double * t, int op, nn_ptr b, nn_srcptr a, slong n,
          const nmod_geometric_progression_struct * G0,
          const nmod_geometric_progression_struct * G1, slong ms)
{
    timeit_t T;
    slong reps = 1;

    for (;;)
    {
        timeit_start_us(T);
        for (slong i = 0; i < reps; i++)
            run_op(op, b, a, n, G0);
        timeit_stop_us(T);
        if (T->wall >= 1000 * ms)
            break;
        reps *= 2;
    }

    t[0] = t[1] = 1e300;
    for (int r = 0; r < 3; r++)
    {
        for (int k = 0; k < 2; k++)
        {
            timeit_start_us(T);
            for (slong i = 0; i < reps; i++)
                run_op(op, b, a, n, k ? G1 : G0);
            timeit_stop_us(T);
            t[k] = FLINT_MIN(t[k], (double) T->wall / reps);
        }
    }
}

/* time in us of one init of length len with function f */
static double
time_init(ulong r, slong len, nmod_t mod, ulong f)
{
    double t;

    TIME_MIN(t, 30, nmod_geometric_progression_t T;
             _nmod_geometric_progression_init_function(T, r, len, mod, f);
             nmod_geometric_progression_clear(T));

    return t;
}

static ulong
modulus(const char * s, flint_rand_t state)
{
    if (!strcmp(s, "fft"))
#if FLINT_BITS == 64
        return UWORD(1108307720798209);  /* 63*2^44 + 1 */
#else
        return UWORD(2013265921);        /* 15*2^27 + 1 */
#endif
    if (!strcmp(s, "ntt"))
        return UWORD(998244353);
    return n_randprime(state, atoi(s), 1);
}

/* same outputs with G0 and G1 */
static void
check(int op, slong n, nn_srcptr a, const nmod_geometric_progression_struct * G0,
      const nmod_geometric_progression_struct * G1)
{
    nn_ptr b = _nmod_vec_init(n), c = _nmod_vec_init(n);
    slong m = (op == 2) ? n / 2 : n;

    run_op(op, b, a, n, G0);
    run_op(op, c, a, n, G1);
    if (!_nmod_vec_equal(b, c, m))
    {
        flint_printf("FAIL: %s, n = %wd\n", names[op], n);
        flint_abort();
    }

    _nmod_vec_clear(b);
    _nmod_vec_clear(c);
}

/* measurements at one length */
typedef struct
{
    double t[3][2];         /* operations, without / with transforms */
    int stored[3];          /* whether the transform of the operation is stored */
    double init[NPRE][2];   /* precomputations, without / with transforms */
}
row_struct;

static void
measure(row_struct * R, slong n, nmod_t mod, ulong r, flint_rand_t state)
{
    nmod_geometric_progression_t G[2];
    nn_ptr a = _nmod_vec_init(n), b = _nmod_vec_init(n);

    _nmod_vec_rand(a, state, n, mod);

    _nmod_geometric_progression_init_function(G[0], r, n, mod, 7);
    _nmod_geometric_progression_init_function(G[1], r, n, mod, 15);

    R->stored[0] = (G[1]->ev_f_fft != NULL);
    R->stored[1] = (G[1]->int_f_fft != NULL);
    R->stored[2] = (G[1]->ext_ff_fft != NULL);

    for (int op = 0; op < 3; op++)
    {
        check(op, n, a, G[0], G[1]);
        time_pair(R->t[op], op, b, a, n, G[0], G[1], 30);
    }

    for (int c = 0; c < NPRE; c++)
        for (int k = 0; k < 2; k++)
            R->init[c][k] = time_init(r, n, mod, pre[c].f | (k ? 8 : 0));

    nmod_geometric_progression_clear(G[0]);
    nmod_geometric_progression_clear(G[1]);
    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
}

/* b/e of precomputation c, see above */
static void
breakeven(char * s, const row_struct * R, int c)
{
    double gain = 0;
    int stored = 0;

    for (int op = 0; op < 3; op++)
    {
        if (((pre[c].f >> op) & 1) && R->stored[op])
        {
            stored = 1;
            gain += R->t[op][0] - R->t[op][1];
        }
    }

    if (!stored)
        strcpy(s, "-");
    else if (gain <= 0)
        strcpy(s, "inf");
    else
        flint_sprintf(s, "%wd", (slong) FLINT_MAX(1.0,
                      ceil((R->init[c][1] - R->init[c][0]) / gain)));
}

static void
table(nmod_t mod, ulong r, int concise, flint_rand_t state)
{
    slong nrows = sizeof(lens) / sizeof(lens[0]), m = 0;
    row_struct * R = flint_malloc(nrows * sizeof(row_struct));
    char be[24];

    /* operations; in concise mode, also the precomputations */
    flint_printf("%6s", "n");
    for (int op = 0; op < 3; op++)
    {
        if (concise)
            flint_printf(" %6s", names[op]);
        else
            flint_printf(" | %8s %8s %5s", names[op], "+tf", "ratio");
    }
    if (concise)
    {
        flint_printf(" |");
        for (int c = 0; c < NPRE; c++)
            flint_printf(" %6s %5s", pre[c].name, "b/e");
    }
    flint_printf("\n");

    for (m = 0; m < nrows; m++)
    {
        slong n = lens[m];

        if (mod.n - 1 < (ulong) 2 * n)
            break;

        measure(R + m, n, mod, r, state);

        flint_printf("%6wd", n);
        for (int op = 0; op < 3; op++)
        {
            double * t = R[m].t[op];

            if (concise)
                flint_printf(" %6.2f", t[0] / t[1]);
            else
                flint_printf(" | %8.2f %8.2f %5.2f", t[0], t[1], t[0] / t[1]);
        }
        if (concise)
        {
            flint_printf(" |");
            for (int c = 0; c < NPRE; c++)
            {
                breakeven(be, R + m, c);
                flint_printf(" %6.2f %5s", R[m].init[c][1] / R[m].init[c][0], be);
            }
        }
        flint_printf("\n");
        fflush(stdout);
    }

    /* precomputations */
    if (!concise)
    {
        flint_printf("\ninit times in us\n%6s", "n");
        for (int c = 0; c < NPRE; c++)
            flint_printf(" | %8s %8s %5s", pre[c].name, "+tf", "b/e");
        flint_printf("\n");

        for (slong i = 0; i < m; i++)
        {
            flint_printf("%6wd", lens[i]);
            for (int c = 0; c < NPRE; c++)
            {
                breakeven(be, R + i, c);
                flint_printf(" | %8.2f %8.2f %5s", R[i].init[c][0], R[i].init[c][1], be);
            }
            flint_printf("\n");
        }
    }

    flint_free(R);
}

/* grid of lengths up to L, ending at L */
static slong
next_len(slong n, slong L)
{
    return (n == L) ? L + 1 : FLINT_MIN(L, n + FLINT_MAX(1, n / 32));
}

static void
sweep(slong L, nmod_t mod, ulong r, flint_rand_t state)
{
    nmod_geometric_progression_t G[2];
    nn_ptr a, b;

    if (mod.n - 1 < (ulong) 2 * L)
    {
        flint_printf("modulus too small for L = %wd\n", L);
        return;
    }

    a = _nmod_vec_init(L);
    b = _nmod_vec_init(L);
    _nmod_vec_rand(a, state, L, mod);
    _nmod_geometric_progression_init_function(G[0], r, L, mod, 7);
    _nmod_geometric_progression_init_function(G[1], r, L, mod, 15);

    for (int op = 0; op < 3; op++)
    {
        double rmin = 1e300, rmax = 0, loss = 0;
        slong nmin = 0, nmax = 0, nloss = 0, count = 0, slow = 0;
        for (slong n = (op == 2) ? 2 : 1; n <= L; n = next_len(n, L))
        {
            double t[2], ratio;

            check(op, n, a, G[0], G[1]);

            time_pair(t, op, b, a, n, G[0], G[1], 10);
            ratio = t[0] / t[1];

            if (ratio < 0.97)
            {
                time_pair(t, op, b, a, n, G[0], G[1], 30);
                ratio = t[0] / t[1];

                if (ratio < 0.97)
                {
                    flint_printf("%6s n = %6wd: %9.2f %9.2f  ratio %5.2f\n",
                                 names[op], n, t[0], t[1], ratio);
                    slow++;
                }
            }

            if (ratio < rmin)
                rmin = ratio, nmin = n;
            if (ratio > rmax)
                rmax = ratio, nmax = n;
            if (t[1] - t[0] > loss)
                loss = t[1] - t[0], nloss = n;
            count++;
        }

        flint_printf("%6s: ratio min %5.2f (n = %wd), max %5.2f (n = %wd), "
                     "< 0.97 at %wd of %wd lengths, largest loss %.3f us (n = %wd)\n",
                     names[op], rmin, nmin, rmax, nmax, slow, count, loss, nloss);
        fflush(stdout);
    }

    nmod_geometric_progression_clear(G[0]);
    nmod_geometric_progression_clear(G[1]);
    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
}

int main(int argc, char * argv[])
{
    const char * defaults[] = {"20", "30", "fft", "ntt", "60", "64"};
    const char ** mods = defaults;
    slong nmods = 6, L = 0;
    int concise = 0, i = 1;
    flint_rand_t state;

    for ( ; i < argc; i++)
    {
        if (!strcmp(argv[i], "--concise"))
            concise = 1;
        else if (!strcmp(argv[i], "-s") && i + 1 < argc)
            L = atol(argv[++i]);
        else
            break;
    }

    if (i < argc)
    {
        mods = (const char **) argv + i;
        nmods = argc - i;
    }

    flint_rand_init(state);

    for (slong m = 0; m < nmods; m++)
    {
        nmod_t mod;

        nmod_init(&mod, modulus(mods[m], state));

        flint_printf("\nmodulus %s (%wu), ratio = without / with transforms",
                     mods[m], mod.n);
        if (L > 0)
            flint_printf(", G of length %wd, times in us", L);
        else if (concise)
            flint_printf(" (init: with / without)");
        else
            flint_printf(", times of one call in us");
        flint_printf("\n");

        if (L > 0)
            sweep(L, mod, n_primitive_root_prime(mod.n), state);
        else
            table(mod, n_primitive_root_prime(mod.n), concise, state);
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
