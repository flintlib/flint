/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Benchmark FLINT's own gemm kernels (flint_sgemm_fallback,
    flint_dgemm_fallback) over a range of shapes and thread counts,
    against the BLAS wrappers (flint_sgemm_blas, flint_dgemm_blas) where
    FLINT is configured with BLAS. Both are called directly rather than
    through flint_sgemm/flint_dgemm, so the comparison does not depend on
    the value of flint_gemm_use_blas.

    To profile the portable backends instead of the native vectors,
    uncomment one of the FLINT_MACHINE_VECTORS_* overrides at the top of
    machine_vectors/gemm.c and rebuild that file.

        p-gemm [options]

          -s          single precision only
          -d          double precision only
          -threads t  comma separated thread counts (default 1,2,4,8)
          -set name   shape set: pow2, round, odd, thin, all (default all)
          -max d      skip shapes with any dimension above d
          -reps r     minimum repetitions per measurement (default 3)

    The BLAS thread count is set to match FLINT's, so that the two are
    compared at equal parallelism. This uses vendor entry points looked
    up as weak symbols (with dlsym on macOS), since there is no portable
    way to do it; if none is found, the program says so, and the
    environment variable for the BLAS in use (e.g. OPENBLAS_NUM_THREADS,
    or VECLIB_MAXIMUM_THREADS for Accelerate) must be set to get a
    meaningful comparison.
*/

#include <string.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "machine_vectors.h"
#include "thread_support.h"
#include "ulong_extras.h"

/*
    Vendor thread controls, declared weak so that this links whether or
    not the BLAS in use provides them (and when there is no BLAS at all).

    This does not work on macOS: there, a weak reference to a function must
    still be resolved by some library at (static) link time, so linking
    fails as soon as one of these entry points is missing from the BLAS in
    use (e.g. mkl_set_num_threads with OpenBLAS, or all of them with
    Accelerate). Look them up at run time with dlsym instead.
*/
#if FLINT_USES_BLAS && defined(__APPLE__)
#include <dlfcn.h>
# define HAVE_DLSYM_BLAS_THREADS 1
#elif defined(__GNUC__) && FLINT_USES_BLAS
extern void openblas_set_num_threads(int) __attribute__((weak));
extern void goto_set_num_threads(int) __attribute__((weak));
extern void mkl_set_num_threads(int) __attribute__((weak));
extern void bli_thread_set_num_threads(long) __attribute__((weak));
# define HAVE_WEAK_BLAS_THREADS 1
#endif

static int
blas_set_num_threads(slong n)
{
#if defined(HAVE_DLSYM_BLAS_THREADS)
    static const char * const int_names[] = {
        "openblas_set_num_threads",
        "goto_set_num_threads",
        "mkl_set_num_threads",
    };
    void (*set_int)(int);
    void (*set_long)(long);
    size_t i;

    for (i = 0; i < sizeof(int_names) / sizeof(int_names[0]); i++)
    {
        *(void **) (&set_int) = dlsym(RTLD_DEFAULT, int_names[i]);
        if (set_int != NULL)
        {
            set_int((int) n);
            return 1;
        }
    }

    *(void **) (&set_long) = dlsym(RTLD_DEFAULT, "bli_thread_set_num_threads");
    if (set_long != NULL)
    {
        set_long((long) n);
        return 1;
    }
#elif defined(HAVE_WEAK_BLAS_THREADS)
    if (openblas_set_num_threads != NULL)
    {
        openblas_set_num_threads((int) n);
        return 1;
    }

    if (goto_set_num_threads != NULL)
    {
        goto_set_num_threads((int) n);
        return 1;
    }

    if (mkl_set_num_threads != NULL)
    {
        mkl_set_num_threads((int) n);
        return 1;
    }

    if (bli_thread_set_num_threads != NULL)
    {
        bli_thread_set_num_threads((long) n);
        return 1;
    }
#else
    (void) n;
#endif

    return 0;
}

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

static double
gflops(slong m, slong k, slong n, double t)
{
    if (t <= 0.0)
        return 0.0;

    return 2.0 * (double) m * (double) k * (double) n / t * 1e-9;
}

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

/* flint_printf does not take width flags with %wd */
static void
print_slong_padded(slong x, slong width)
{
    slong digits = 1, tmp = x;

    while (tmp >= 10)
    {
        tmp /= 10;
        digits++;
    }

    flint_printf("%wd", x);

    while (digits++ < width)
        flint_printf(" ");
}

int
main(int argc, char ** argv)
{
    slong i, j;
    int do_single = 1, do_double = 1;
    slong nthreads[16];
    slong num_thread_counts = 4;
    const char * set = "all";
    slong maxdim = WORD_MAX;
    slong reps = 3;
    flint_rand_t state;

    nthreads[0] = 1; nthreads[1] = 2; nthreads[2] = 4; nthreads[3] = 8;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-s"))
            do_double = 0;
        else if (!strcmp(argv[i], "-d"))
            do_single = 0;
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
            flint_printf("usage: %s [-s] [-d] [-threads 1,2,4,8] "
                         "[-set pow2|round|odd|thin|all] [-max d] "
                         "[-reps r]\n", argv[0]);
            return 1;
        }
    }

    flint_rand_init(state);

#if FLINT_USES_BLAS
    flint_printf("flint gemm kernels vs BLAS (FLINT_USES_BLAS)\n");
    flint_printf("ratio = blas time / flint time; > 1 means flint is faster\n");

    if (!blas_set_num_threads(1))
        flint_printf("WARNING: cannot set the BLAS thread count from here; "
                     "set it in the environment (e.g. OPENBLAS_NUM_THREADS) "
                     "to match, or the comparison is meaningless\n");
#else
    flint_printf("no BLAS configured; timing the flint kernels only\n");
#endif
    flint_printf("times are best of >= %wd runs, in ms\n\n", reps);

    for (i = 0; i < (slong) NUM_SHAPES; i++)
    {
        slong m = shapes[i].m, k = shapes[i].k, n = shapes[i].n;
        float *Af, *Bf, *Cf;
        double *Ad, *Bd, *Cd;
        slong e;

        if (strcmp(set, "all") != 0 && strcmp(set, shapes[i].set) != 0)
            continue;

        if (m > maxdim || k > maxdim || n > maxdim)
            continue;

        Af = flint_malloc(m * k * sizeof(float));
        Bf = flint_malloc(k * n * sizeof(float));
        Cf = flint_malloc(m * n * sizeof(float));
        Ad = flint_malloc(m * k * sizeof(double));
        Bd = flint_malloc(k * n * sizeof(double));
        Cd = flint_malloc(k * n * sizeof(double));
        Cd = flint_realloc(Cd, m * n * sizeof(double));

        for (e = 0; e < m * k; e++)
        {
            Ad[e] = (double) n_randint(state, 100) / 8.0;
            Af[e] = (float) Ad[e];
        }

        for (e = 0; e < k * n; e++)
        {
            Bd[e] = (double) n_randint(state, 100) / 8.0;
            Bf[e] = (float) Bd[e];
        }

        flint_printf("%5wd x %5wd x %5wd  [%s]\n", m, k, n, shapes[i].set);
        flint_printf("  %-8s", "threads");
        if (do_single)
            flint_printf("  %10s %9s", "sgemm", "GFLOP/s");
#if FLINT_USES_BLAS
        if (do_single)
            flint_printf("  %10s %9s %7s", "blas_s", "GFLOP/s", "ratio");
#endif
        if (do_double)
            flint_printf("  %10s %9s", "dgemm", "GFLOP/s");
#if FLINT_USES_BLAS
        if (do_double)
            flint_printf("  %10s %9s %7s", "blas_d", "GFLOP/s", "ratio");
#endif
        flint_printf("\n");

        for (j = 0; j < num_thread_counts; j++)
        {
            slong t = nthreads[j];
            double tf = 0.0, td = 0.0, tbf = 0.0, tbd = 0.0;

            flint_set_num_threads(t);
#if FLINT_USES_BLAS
            blas_set_num_threads(t);
#endif

            if (do_single)
                TIME_BEST(tf, reps,
                    flint_sgemm_fallback(m, k, n, Af, k, Bf, n, Cf, n));

            if (do_double)
                TIME_BEST(td, reps,
                    flint_dgemm_fallback(m, k, n, Ad, k, Bd, n, Cd, n));

#if FLINT_USES_BLAS
            if (do_single)
                TIME_BEST(tbf, reps,
                    flint_sgemm_blas(m, k, n, Af, k, Bf, n, Cf, n));

            if (do_double)
                TIME_BEST(tbd, reps,
                    flint_dgemm_blas(m, k, n, Ad, k, Bd, n, Cd, n));
#endif

            flint_printf("  ");
            print_slong_padded(t, 8);

            if (do_single)
                flint_printf("  %10.4f %9.2f", tf * 1e3, gflops(m, k, n, tf));
#if FLINT_USES_BLAS
            if (do_single)
                flint_printf("  %10.4f %9.2f %7.2f", tbf * 1e3,
                             gflops(m, k, n, tbf),
                             tbf > 0.0 ? tbf / tf : 0.0);
#endif
            if (do_double)
                flint_printf("  %10.4f %9.2f", td * 1e3, gflops(m, k, n, td));
#if FLINT_USES_BLAS
            if (do_double)
                flint_printf("  %10.4f %9.2f %7.2f", tbd * 1e3,
                             gflops(m, k, n, tbd),
                             tbd > 0.0 ? tbd / td : 0.0);
#endif
            flint_printf("\n");
        }

        flint_printf("\n");

        flint_free(Af); flint_free(Bf); flint_free(Cf);
        flint_free(Ad); flint_free(Bd); flint_free(Cd);
    }

    flint_set_num_threads(1);
    flint_rand_clear(state);
    flint_cleanup_master();

    return 0;
}
