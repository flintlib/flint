/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5 and Opus 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Benchmark the small-modulus matrix multiplication over a range of
    moduli, shapes and thread counts:

      u8    nmod_mat_mul_u8, the public entry, including the
            ulong <-> byte conversions;
      raw   _nmod_mat_mul_u8 on byte matrices, the core without them;
      blas  nmod_mat_mul_blas, the baseline, backed by FLINT's own
            flint_sgemm when no external BLAS is configured.

    The shapes are those of machine_vectors/profile/p-gemm.c.

        p-mul_u8 [options]

          -p x,y,...  comma separated moduli (default 2,3,5,13,127,251)
          -fn a,b,..  subset of u8,raw,blas,classical (default the
                      first three; classical = nmod_mat_mul_classical,
                      the reference for small-size crossovers)
          -threads t  comma separated thread counts (default 1,2,4,8)
          -set name   shape set: pow2, round, odd, thin, all (default all)
          -max d      skip shapes with any dimension above d
          -reps r     minimum repetitions per measurement (default 3)

    Each (modulus, shape) first cross-checks the selected functions
    against nmod_mat_mul_classical on the same operands, so a timing
    table cannot be produced from wrong results. With an external BLAS
    configured, its thread count is not under flint_set_num_threads;
    a warning is printed in that case.
*/

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "profiler.h"
#include "flint.h"
#include "nmod.h"
#include "nmod_mat.h"
#include "thread_support.h"
#include "ulong_extras.h"

typedef struct
{
    slong m, k, n;
    const char * set;
}
shape_t;

static const shape_t shapes[] =
{
    /* all powers of two up to 4096 */
    {    8,    8,    8, "pow2" },
    {   16,   16,   16, "pow2" },
    {   32,   32,   32, "pow2" },
    {   64,   64,   64, "pow2" },
    {  128,  128,  128, "pow2" },
    {  256,  256,  256, "pow2" },
    {  512,  512,  512, "pow2" },
    { 1024, 1024, 1024, "pow2" },
    { 2048, 2048, 2048, "pow2" },
    { 4096, 4096, 4096, "pow2" },

    /* multiples of 10, 100, 1000 */
    {   10,   10,   10, "round" },
    {   50,   50,   50, "round" },
    {  100,  100,  100, "round" },
    {  300,  300,  300, "round" },
    {  500,  500,  500, "round" },
    { 1000, 1000, 1000, "round" },
    { 2000, 2000, 2000, "round" },
    { 3000, 3000, 3000, "round" },

    /* odd and near-power-of-two sizes */
    {    7,    7,    7, "odd" },
    {   17,   17,   17, "odd" },
    {   31,   33,   29, "odd" },
    {  101,  103,   97, "odd" },
    {  255,  257,  255, "odd" },
    {  999, 1001, 1003, "odd" },
    { 1234,  777, 1919, "odd" },
    { 2047, 2049, 2048, "odd" },

    /* unbalanced shapes */
    {    1, 2048, 2048, "thin" },
    { 2048, 2048,    1, "thin" },
    {    8, 2048, 2048, "thin" },
    { 2048, 2048,    8, "thin" },
    { 2048,    8, 2048, "thin" },
    {   64, 4096, 4096, "thin" },
    { 4096,   64, 4096, "thin" },
    { 4096, 4096,   64, "thin" },
    {  100, 5000,  300, "thin" },
};

#define NUM_SHAPES (sizeof(shapes) / sizeof(shape_t))

/*
    Best of several timing windows, each running enough repetitions to
    exceed 50 ms; microsecond timers, since a single small multiplication
    is far below millisecond resolution.
*/
#define TIME_BEST(t, minreps, expr) \
    do { \
        slong __reps = (minreps), __i, __w; \
        double __best = 1e300, __cur; \
        timeit_t __timer; \
        for (;;) \
        { \
            timeit_start_us(__timer); \
            for (__i = 0; __i < __reps; __i++) { expr; } \
            timeit_stop_us(__timer); \
            if (__timer->wall >= 50000 || __reps >= WORD(1) << 22) \
                break; \
            __reps *= 4; \
        } \
        __best = (double) __timer->wall * 1e-6 / (double) __reps; \
        for (__w = 0; __w < 2; __w++) \
        { \
            timeit_start_us(__timer); \
            for (__i = 0; __i < __reps; __i++) { expr; } \
            timeit_stop_us(__timer); \
            __cur = (double) __timer->wall * 1e-6 / (double) __reps; \
            if (__cur < __best) \
                __best = __cur; \
        } \
        (t) = __best; \
    } while (0)

#define FN_U8        0
#define FN_RAW       1
#define FN_BLAS      2
#define FN_CLASSICAL 3

static const char * fn_names[4] = { "u8", "raw", "blas", "classical" };

int
main(int argc, char ** argv)
{
    slong i, j, t, f;
    ulong moduli[16] = { 2, 3, 5, 13, 127, 251 };
    slong num_moduli = 6;
    int do_fn[4] = { 1, 1, 1, 0 };
    slong nthreads[16];
    slong num_thread_counts = 4;
    const char * set = "all";
    slong maxdim = WORD_MAX;
    slong reps = 3;
    flint_rand_t state;

    nthreads[0] = 1; nthreads[1] = 2; nthreads[2] = 4; nthreads[3] = 8;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-p") && i + 1 < argc)
        {
            char * s = argv[++i];

            num_moduli = 0;
            while (*s != '\0' && num_moduli < 16)
            {
                moduli[num_moduli++] = (ulong) atol(s);
                while (*s != '\0' && *s != ',')
                    s++;
                if (*s == ',')
                    s++;
            }
        }
        else if (!strcmp(argv[i], "-fn") && i + 1 < argc)
        {
            char * s = argv[++i];

            do_fn[0] = do_fn[1] = do_fn[2] = do_fn[3] = 0;
            while (*s != '\0')
            {
                if (!strncmp(s, "u8", 2) && (s[2] == ',' || s[2] == '\0'))
                    do_fn[FN_U8] = 1;
                else if (!strncmp(s, "raw", 3))
                    do_fn[FN_RAW] = 1;
                else if (!strncmp(s, "blas", 4))
                    do_fn[FN_BLAS] = 1;
                else if (!strncmp(s, "classical", 9))
                    do_fn[FN_CLASSICAL] = 1;
                while (*s != '\0' && *s != ',')
                    s++;
                if (*s == ',')
                    s++;
            }
        }
        else if (!strcmp(argv[i], "-threads") && i + 1 < argc)
        {
            char * s = argv[++i];

            num_thread_counts = 0;
            while (*s != '\0' && num_thread_counts < 16)
            {
                nthreads[num_thread_counts++] = atol(s);
                while (*s != '\0' && *s != ',')
                    s++;
                if (*s == ',')
                    s++;
            }
        }
        else if (!strcmp(argv[i], "-set") && i + 1 < argc)
            set = argv[++i];
        else if (!strcmp(argv[i], "-max") && i + 1 < argc)
            maxdim = atol(argv[++i]);
        else if (!strcmp(argv[i], "-reps") && i + 1 < argc)
            reps = atol(argv[++i]);
        else
        {
            flint_printf("usage: %s [-p 2,3,5,13,127,251] [-fn u8,raw,blas,classical] "
                         "[-threads 1,2,4,8] [-set pow2|round|odd|thin|all] "
                         "[-max d] [-reps r]\n", argv[0]);
            return 1;
        }
    }

#if FLINT_USES_BLAS
    if (do_fn[FN_BLAS])
        flint_printf("WARNING: FLINT is configured to use an external "
                     "BLAS which may not respect flint_set_num_threads. "
                     "Set the BLAS thread count in the environment (e.g. "
                     "OPENBLAS_NUM_THREADS) to match, or the comparison "
                     "is meaningless.\n\n");
#endif

    flint_rand_init(state);

    flint_printf("u8 = nmod_mat_mul_u8, raw = _nmod_mat_mul_u8 (byte "
                 "matrices, no ulong conversion), blas = nmod_mat_mul_blas\n");
    flint_printf("times are best of >= %wd runs, in ms; - means the "
                 "function declined the problem (blas only)\n\n", reps);

    /* header */
    flint_printf("%-20s %-4s %-5s", "m x k x n", "p", "fn");
    for (t = 0; t < num_thread_counts; t++)
        flint_printf("   nt=%wd     ", nthreads[t]);
    flint_printf("\n");

    for (i = 0; i < (slong) NUM_SHAPES; i++)
    {
        slong m = shapes[i].m, k = shapes[i].k, n = shapes[i].n;
        nmod_mat_t A, B, C, D;
        uint8_t *Au, *Bu, *Cu, *Du;
        slong e;

        if (strcmp(set, "all") != 0 && strcmp(set, shapes[i].set) != 0)
            continue;

        if (m > maxdim || k > maxdim || n > maxdim)
            continue;

        Au = flint_malloc(m * k);
        Bu = flint_malloc(k * n);
        Cu = flint_malloc(m * n);
        Du = flint_malloc(m * n);

        for (j = 0; j < num_moduli; j++)
        {
            ulong p = moduli[j];
            int fn_ok[4] = { 0, 0, 0, 1 };
            nmod_t mod;

            nmod_init(&mod, p);

            nmod_mat_init(A, m, k, p);
            nmod_mat_init(B, k, n, p);
            nmod_mat_init(C, m, n, p);
            nmod_mat_init(D, m, n, p);

            nmod_mat_randfull(A, state);
            nmod_mat_randfull(B, state);

            for (e = 0; e < m * k; e++)
                Au[e] = (uint8_t) nmod_mat_entry(A, e / k, e % k);
            for (e = 0; e < k * n; e++)
                Bu[e] = (uint8_t) nmod_mat_entry(B, e / n, e % n);

            /*
                Cross-check every selected function on these operands at
                the largest requested thread count before timing anything.
            */
            flint_set_num_threads(nthreads[num_thread_counts - 1]);
            nmod_mat_mul_classical(D, A, B);

            if (do_fn[FN_U8])
            {
                nmod_mat_zero(C);
                nmod_mat_mul_u8(C, A, B);
                fn_ok[FN_U8] = 1;
                if (!nmod_mat_equal(C, D))
                {
                    flint_printf("FAIL: nmod_mat_mul_u8 disagrees with "
                                 "classical at p = %wu, %wd x %wd x %wd\n",
                                 p, m, k, n);
                    flint_abort();
                }
            }

            if (do_fn[FN_RAW])
            {
                memset(Cu, 0xee, m * n);
                _nmod_mat_mul_u8(Cu, n, Au, k, Bu, n, m, k, n, mod);
                fn_ok[FN_RAW] = 1;
                {
                    for (e = 0; e < m * n; e++)
                        Du[e] = (uint8_t) nmod_mat_entry(D, e / n, e % n);
                    if (memcmp(Cu, Du, m * n) != 0)
                    {
                        flint_printf("FAIL: _nmod_mat_mul_u8 disagrees with "
                                     "classical at p = %wu, %wd x %wd x %wd\n",
                                     p, m, k, n);
                        flint_abort();
                    }
                }
            }

            if (do_fn[FN_BLAS])
            {
                nmod_mat_zero(C);
                fn_ok[FN_BLAS] = nmod_mat_mul_blas(C, A, B);
                if (fn_ok[FN_BLAS] && !nmod_mat_equal(C, D))
                {
                    flint_printf("FAIL: nmod_mat_mul_blas disagrees with "
                                 "classical at p = %wu, %wd x %wd x %wd\n",
                                 p, m, k, n);
                    flint_abort();
                }
            }

            for (f = 0; f < 4; f++)
            {
                double tm;

                if (!do_fn[f])
                    continue;

                flint_printf("%wd x %wd x %wd", m, k, n);
                {
                    /* pad the shape column to 20 columns by hand:
                       flint_printf has no %*d */
                    slong len = 0, tmp;
                    for (tmp = m; tmp > 0 || len == 0; tmp /= 10) len++;
                    for (tmp = k; tmp > 0; tmp /= 10) len++;
                    for (tmp = n; tmp > 0; tmp /= 10) len++;
                    len += 6;
                    while (len++ < 20)
                        flint_printf(" ");
                }
                flint_printf(" %-4wu %-5s", p, fn_names[f]);

                for (t = 0; t < num_thread_counts; t++)
                {
                    if (!fn_ok[f])
                    {
                        flint_printf("     -       ");
                        continue;
                    }

                    flint_set_num_threads(nthreads[t]);

                    if (f == FN_U8)
                        TIME_BEST(tm, reps, nmod_mat_mul_u8(C, A, B));
                    else if (f == FN_RAW)
                        TIME_BEST(tm, reps, _nmod_mat_mul_u8(Cu, n, Au, k,
                                                             Bu, n, m, k, n,
                                                             mod));
                    else if (f == FN_BLAS)
                        TIME_BEST(tm, reps, nmod_mat_mul_blas(C, A, B));
                    else
                        TIME_BEST(tm, reps, nmod_mat_mul_classical(C, A, B));

                    flint_printf("%10.4f  ", tm * 1e3);
                }

                flint_printf("\n");
            }

            nmod_mat_clear(A);
            nmod_mat_clear(B);
            nmod_mat_clear(C);
            nmod_mat_clear(D);
        }

        flint_free(Au);
        flint_free(Bu);
        flint_free(Cu);
        flint_free(Du);

        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
