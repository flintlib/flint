/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Benchmark the SIMD matrix multiplications nmod_mat_mul_u32 (moduli
    below 2^32), its uint32-entry variant _nmod_mat_mul_u32,
    nmod_mat_mul_u52 (AVX512-IFMA, moduli up to 2^52), nmod_mat_mul_k52
    (two-limb integer Karatsuba, any SIMD, moduli up to 2^52) and
    nmod_mat_mul_fp50 (double precision with the mulmod of fft_small,
    moduli below 2^50) against the other matrix multiplications of
    nmod_mat, over a range of modulus sizes, dimensions and thread counts,
    and report where a SIMD kernel is the fastest. The purpose is to tune
    the dispatch in nmod_mat_mul, whose parameters live in flint-mparam.h
    (FLINT_NMOD_MAT_MUL_U32_*, _U52_*, _K52_* and _FP50_*).

        p-mul_tune [options]

          -bits b,...  modulus bit sizes; for each, the modulus is a
                       fixed "generic" prime of that size (with a large
                       2^32 mod n, i.e. the slower two-round folding);
                       default 16,20,22,24,26,28,29,30,31,32,36,40,44,
                       48,50,52
          -p n,...     explicit moduli instead (any n up to 2^52)
          -fn f,...    subset of u32,nmod32,u52,k52,fp50,blas,fgemm,
                       classical,strassen,threaded,mul (default: all
                       available;
                       nmod32 is _nmod_mat_mul_u32 on uint32 arrays
                       converted beforehand, threaded is
                       nmod_mat_mul_classical_threaded, mul the dispatch
                       nmod_mat_mul, and fgemm is nmod_mat_mul_blas forced
                       onto FLINT's own gemm, which only differs from blas
                       when FLINT is built with an external BLAS)
          -dims d,...  square dimensions (default 32,40,48,56,64,80,96,
                       112,128,160,192,224,256,320,384,448,512,768,1024,
                       1536,2048,4096)
          -shape m,k,n a single rectangular shape instead of -dims
          -threads t,..thread counts (default 1)
          -tmax s      once a function takes more than s seconds on some
                       dimension it is skipped for the larger ones
                       (default 4)
          -reps r      minimum repetitions per measurement (default 3)
          -csv         additionally print one machine-readable line per
                       measurement: "csv,bits,modulus,threads,m,k,n,fn,us"

    Every (modulus, shape) first cross-checks the SIMD kernels against
    nmod_mat_mul_classical (for dimensions up to 600) so that a timing
    table cannot come from wrong results. Each table row ends with the
    name of the fastest function overall and the ratio
    time(best SIMD kernel)/time(best base algorithm), where the SIMD
    kernels are u32, nmod32, u52, k52 and fp50 (those that are enabled and
    handle the modulus) and the base algorithms are blas, fgemm, classical
    and threaded: strassen and mul are excluded from the ratio because they
    recurse through nmod_mat_mul and so already contain whatever the
    dispatch picks as their leaf; comparing against them would say nothing
    about where the kernels themselves should be used. A ratio below 1
    means that some SIMD kernel beats every base algorithm, and the "best"
    column says which. The summary per modulus and thread count gives the
    smallest dimension from which that holds for all larger dimensions
    tested, which is the crossover to put in nmod_mat_mul; the strassen
    and mul columns then show whether a Strassen level on top pays and
    whether the dispatch currently makes the right choice.

    The dimensions are dense up to 512 so that the crossovers can be
    read off; 4096 shows the asymptotic regime.

    External BLAS: an OpenMP or pthreads BLAS decides its own thread
    count and ignores flint_set_num_threads, which makes a comparison
    at a given thread count meaningless. This program therefore pins
    the usual variables (OMP_NUM_THREADS, BLIS_NUM_THREADS,
    OPENBLAS_NUM_THREADS, MKL_NUM_THREADS) to 1 at startup unless they
    are already set in the environment. This works for libraries that
    read them at first use (BLIS does); for one that reads them when
    it is loaded, export the variables before launching. With -threads
    beyond 1, the blas column of an external BLAS is thus single
    threaded; fgemm and everything else follow flint_set_num_threads.
*/

/* for setenv under -std=c11 */
#if !defined(_WIN32) && !defined(_MSC_VER) && !defined(_POSIX_C_SOURCE)
# define _POSIX_C_SOURCE 200112L
#endif

#include <string.h>
#include <stdlib.h>
#include "profiler.h"
#include "machine_vectors.h"
#include "flint.h"
#include "nmod.h"
#include "nmod_mat.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod_mat/impl.h"

#define FN_U32       0
#define FN_NMOD32    1
#define FN_U52       2
#define FN_K52       3
#define FN_FP50      4
#define FN_BLAS      5
#define FN_FGEMM     6
#define FN_CLASSICAL 7
#define FN_STRASSEN  8
#define FN_THREADED  9
#define FN_MUL       10
#define NUM_FN       11

static const char * fn_names[NUM_FN] =
    { "u32", "nmod32", "u52", "k52", "fp50", "blas", "fgemm", "classical",
      "strassen", "threaded", "mul" };

/* the SIMD kernels being tuned, and the functions they are compared
   against for the crossover: those that do not themselves go through
   nmod_mat_mul and are not SIMD kernels */
static const int fn_is_simd[NUM_FN] = { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 };
static const int fn_is_base[NUM_FN] = { 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0 };

/* the nmod_mat-level SIMD kernels, for the cross-check and the timing */
typedef int (* simd_fn)(nmod_mat_t, const nmod_mat_t, const nmod_mat_t);
static const simd_fn simd_funcs[NUM_FN] =
    { nmod_mat_mul_u32, NULL, nmod_mat_mul_u52, nmod_mat_mul_k52,
      nmod_mat_mul_fp50, NULL, NULL, NULL, NULL, NULL, NULL };

/* uint32 copies of A, B and C for _nmod_mat_mul_u32 */
typedef struct
{
    uint32_t * a, * b, * c;
    slong m, k, n;
    nmod_t mod;
}
nmod32_data;

static void
nmod32_data_init(nmod32_data * w, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j;

    w->m = A->r; w->k = A->c; w->n = B->c;
    w->mod = A->mod;
    w->a = flint_malloc(FLINT_MAX(w->m * w->k, 1) * sizeof(uint32_t));
    w->b = flint_malloc(FLINT_MAX(w->k * w->n, 1) * sizeof(uint32_t));
    w->c = flint_malloc(FLINT_MAX(w->m * w->n, 1) * sizeof(uint32_t));

    for (i = 0; i < w->m; i++)
        for (j = 0; j < w->k; j++)
            w->a[i * w->k + j] = (uint32_t) nmod_mat_entry(A, i, j);
    for (i = 0; i < w->k; i++)
        for (j = 0; j < w->n; j++)
            w->b[i * w->n + j] = (uint32_t) nmod_mat_entry(B, i, j);
}

static void
nmod32_data_clear(nmod32_data * w)
{
    flint_free(w->a);
    flint_free(w->b);
    flint_free(w->c);
}

static int
nmod32_mul(nmod32_data * w)
{
    return _nmod_mat_mul_u32(w->c, w->n, w->a, w->k, w->b, w->n,
                             w->m, w->k, w->n, w->mod);
}

/* does the uint32 result agree with the nmod_mat D ? */
static int
nmod32_equal(const nmod32_data * w, const nmod_mat_t D)
{
    slong i, j;

    for (i = 0; i < w->m; i++)
        for (j = 0; j < w->n; j++)
            if (w->c[i * w->n + j] != nmod_mat_entry(D, i, j))
                return 0;

    return 1;
}

/* nmod_mat_mul_blas on FLINT's own gemm whatever the build */
static int
nmod_mat_mul_fgemm(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    int save = flint_gemm_use_blas, ret;

    flint_gemm_use_blas = 0;
    ret = nmod_mat_mul_blas(C, A, B);
    flint_gemm_use_blas = save;

    return ret;
}

/* best of several timing windows, each long enough to be measurable */
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
            if (__best * __reps > 2.0) break; \
            timeit_start_us(__timer); \
            for (__i = 0; __i < __reps; __i++) { expr; } \
            timeit_stop_us(__timer); \
            __cur = (double) __timer->wall * 1e-6 / (double) __reps; \
            if (__cur < __best) \
                __best = __cur; \
        } \
        (t) = __best; \
    } while (0)

static slong
parse_list(const char * s, slong * out, slong max)
{
    slong num = 0;

    while (*s != '\0' && num < max)
    {
        out[num++] = atol(s);
        while (*s != '\0' && *s != ',')
            s++;
        if (*s == ',')
            s++;
    }

    return num;
}

/* a prime of the given bit size, away from powers of two so that
   2^32 mod n is large (the harder case for the in-kernel folding) */
static ulong
generic_prime(slong bits)
{
    ulong lo, n;

    if (bits <= 2)
        return (bits <= 1) ? 2 : 3;

    lo = UWORD(1) << (bits - 1);
    n = lo + lo / 2 + lo / 8 + 12345 % lo;   /* about 1.6 * 2^(bits-1) */
    n = n_nextprime(n, 1);

    /* n_nextprime may step past 2^bits for tiny sizes */
    while (bits < FLINT_BITS && n >= (UWORD(1) << bits))
        n = n_nextprime(lo + lo / 8, 1);

    return n;
}

int
main(int argc, char ** argv)
{
    slong i, j, t, f, d;
    slong bits_list[64] = { 16, 20, 22, 24, 26, 28, 29, 30, 31, 32,
                            36, 40, 44, 48, 50, 52 };
    slong num_bits = 16;
    ulong moduli[64];
    slong num_moduli = 0;
    int explicit_moduli = 0;
    /* fgemm would just repeat blas without an external BLAS */
#if FLINT_USES_BLAS
    int do_fn[NUM_FN] = { 1, 1, NMOD_MAT_HAVE_MUL_U52, 1, 1, 1, 1, 1, 1, 1, 1 };
#else
    int do_fn[NUM_FN] = { 1, 1, NMOD_MAT_HAVE_MUL_U52, 1, 1, 1, 0, 1, 1, 1, 1 };
#endif
    slong nthreads[16] = { 1 };
    slong num_thread_counts = 1;
    slong dims[64] = { 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224,
                       256, 320, 384, 448, 512, 768, 1024, 1536, 2048, 4096 };
    slong num_dims = 22;
    slong shape[3] = { 0, 0, 0 };
    int have_shape = 0;
    double tmax = 4.0;
    slong reps = 3;
    int csv = 0;
    flint_rand_t state;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-bits") && i + 1 < argc)
            num_bits = parse_list(argv[++i], bits_list, 64);
        else if (!strcmp(argv[i], "-p") && i + 1 < argc)
        {
            slong tmp[64];
            num_moduli = parse_list(argv[++i], tmp, 64);
            for (j = 0; j < num_moduli; j++)
                moduli[j] = (ulong) tmp[j];
            explicit_moduli = 1;
        }
        else if (!strcmp(argv[i], "-fn") && i + 1 < argc)
        {
            char * s = argv[++i];

            for (f = 0; f < NUM_FN; f++)
                do_fn[f] = 0;

            while (*s != '\0')
            {
                for (f = 0; f < NUM_FN; f++)
                {
                    size_t len = strlen(fn_names[f]);
                    if (!strncmp(s, fn_names[f], len)
                            && (s[len] == ',' || s[len] == '\0'))
                        do_fn[f] = 1;
                }
                while (*s != '\0' && *s != ',')
                    s++;
                if (*s == ',')
                    s++;
            }
        }
        else if (!strcmp(argv[i], "-dims") && i + 1 < argc)
            num_dims = parse_list(argv[++i], dims, 64);
        else if (!strcmp(argv[i], "-shape") && i + 1 < argc)
        {
            if (parse_list(argv[++i], shape, 3) != 3)
            {
                flint_printf("-shape wants m,k,n\n");
                return 1;
            }
            have_shape = 1;
        }
        else if (!strcmp(argv[i], "-threads") && i + 1 < argc)
            num_thread_counts = parse_list(argv[++i], nthreads, 16);
        else if (!strcmp(argv[i], "-tmax") && i + 1 < argc)
            tmax = atof(argv[++i]);
        else if (!strcmp(argv[i], "-reps") && i + 1 < argc)
            reps = atol(argv[++i]);
        else if (!strcmp(argv[i], "-csv"))
            csv = 1;
        else
        {
            flint_printf("usage: %s [-bits 16,20,...,32] [-p n,...] "
                         "[-fn u32,nmod32,u52,k52,fp50,blas,fgemm,classical,strassen,threaded,mul] "
                         "[-dims 32,...,4096] [-shape m,k,n] "
                         "[-threads 1,2,4] [-tmax s] [-reps r] [-csv]\n",
                         argv[0]);
            return 1;
        }
    }

    if (!explicit_moduli)
    {
        num_moduli = num_bits;
        for (j = 0; j < num_bits; j++)
            moduli[j] = generic_prime(bits_list[j]);
    }

    if (have_shape)
        num_dims = 1;

#if FLINT_USES_BLAS
    /* pin an external BLAS to one thread unless the user chose otherwise;
       see the header comment */
    {
        const char * vars[] = { "OMP_NUM_THREADS", "BLIS_NUM_THREADS",
                                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS" };
        int pinned = 0;

        for (i = 0; i < 4; i++)
            if (getenv(vars[i]) == NULL)
            {
#if defined(_WIN32) || defined(_MSC_VER)
                _putenv_s(vars[i], "1");
#else
                setenv(vars[i], "1", 0);
#endif
                pinned = 1;
            }

        flint_printf("FLINT is built with an external BLAS: %s pinned to 1 "
                     "thread via the *_NUM_THREADS variables%s. Its blas "
                     "column ignores -threads; fgemm is the same lifting "
                     "and CRT on FLINT's own gemm, which does follow it.\n\n",
                     pinned ? "it was" : "they were not changed,",
                     pinned ? "" : " already set in the environment");
    }
#endif

    flint_rand_init(state);

    flint_printf("u32 = nmod_mat_mul_u32, nmod32 = _nmod_mat_mul_u32 on uint32 "
                 "arrays, u52 = nmod_mat_mul_u52 (AVX512-IFMA%s),\n"
                 "k52 = nmod_mat_mul_k52 (two-limb integer Karatsuba, moduli up to 2^52), "
                 "fp50 = nmod_mat_mul_fp50 (double precision, moduli below 2^50),\n"
                 "blas = nmod_mat_mul_blas, fgemm = blas on FLINT's own gemm, "
                 "classical = nmod_mat_mul_classical,\n"
                 "strassen = nmod_mat_mul_strassen, threaded = "
                 "nmod_mat_mul_classical_threaded, mul = nmod_mat_mul (the dispatch)\n",
                 NMOD_MAT_HAVE_MUL_U52 ? "" : ", not available in this build");
    flint_printf("times in microseconds, best of >= %wd runs; - means declined or "
                 "skipped (slower than %.1f s at a smaller size);\n"
                 "ratio = time(fastest of u32, nmod32, u52, k52, fp50) / time(fastest of blas, fgemm, classical, threaded);\n"
                 "strassen and mul recurse through nmod_mat_mul and are left out of it\n\n",
                 reps, tmax);

    for (j = 0; j < num_moduli; j++)
    {
        ulong p = moduli[j];
        slong pbits = FLINT_BIT_COUNT(p);

        for (t = 0; t < num_thread_counts; t++)
        {
            int skipped[NUM_FN] = { 0 };
            slong first_win = -1;   /* smallest dim from which a SIMD kernel wins */
            int simd_ever_lost_after = 0;
            slong last_dim_measured = 0;

            flint_set_num_threads(nthreads[t]);

            flint_printf("modulus %wu (%wd bits), %wd thread%s\n",
                         p, pbits, nthreads[t], nthreads[t] > 1 ? "s" : "");
            flint_printf("%-18s", "m x k x n");
            for (f = 0; f < NUM_FN; f++)
                if (do_fn[f])
                    flint_printf(" %10s", fn_names[f]);
            flint_printf("   %-10s %s\n", "best", "ratio");

            for (d = 0; d < num_dims; d++)
            {
                slong m, k, n;
                nmod_mat_t A, B, C, D;
                double tm[NUM_FN];
                int ok[NUM_FN];
                double best_other = 1e300, best_simd = 1e300;
                slong best_fn = -1;
                slong shape_len;
                nmod32_data w32 = { NULL, NULL, NULL, 0, 0, 0, { 0, 0, 0 } };
                int do_nmod32 = 1;

                if (have_shape)
                {
                    m = shape[0]; k = shape[1]; n = shape[2];
                }
                else
                    m = k = n = dims[d];

                nmod_mat_init(A, m, k, p);
                nmod_mat_init(B, k, n, p);
                nmod_mat_init(C, m, n, p);
                nmod_mat_randfull(A, state);
                nmod_mat_randfull(B, state);

                if (do_fn[FN_NMOD32] && pbits <= 32)
                    nmod32_data_init(&w32, A, B);
                else
                    do_nmod32 = 0;

                /* cross-check the SIMD kernels while classical is affordable */
                if (FLINT_MAX(FLINT_MAX(m, k), n) <= 600)
                {
                    nmod_mat_init(D, m, n, p);
                    nmod_mat_mul_classical(D, A, B);
                    for (f = 0; f < NUM_FN; f++)
                    {
                        if (!do_fn[f] || simd_funcs[f] == NULL)
                            continue;
                        nmod_mat_randtest(C, state);
                        if (simd_funcs[f](C, A, B) && !nmod_mat_equal(C, D))
                        {
                            flint_printf("FAIL: %s disagrees with classical "
                                         "at p = %wu, %wd x %wd x %wd, "
                                         "%wd threads\n", fn_names[f], p,
                                         m, k, n, nthreads[t]);
                            flint_abort();
                        }
                    }
                    if (do_nmod32 && nmod32_mul(&w32) && !nmod32_equal(&w32, D))
                    {
                        flint_printf("FAIL: _nmod_mat_mul_u32 disagrees with "
                                     "classical at p = %wu, %wd x %wd x %wd, "
                                     "%wd threads\n", p, m, k, n, nthreads[t]);
                        flint_abort();
                    }
                    nmod_mat_clear(D);
                }

                for (f = 0; f < NUM_FN; f++)
                {
                    ok[f] = 0;
                    tm[f] = 0.0;

                    if (!do_fn[f] || skipped[f])
                        continue;

                    switch (f)
                    {
                        case FN_U32:
                            ok[f] = nmod_mat_mul_u32(C, A, B);
                            if (ok[f])
                                TIME_BEST(tm[f], reps, nmod_mat_mul_u32(C, A, B));
                            break;
                        case FN_NMOD32:
                            ok[f] = do_nmod32 && nmod32_mul(&w32);
                            if (ok[f])
                                TIME_BEST(tm[f], reps, nmod32_mul(&w32));
                            break;
                        case FN_U52:
                            ok[f] = nmod_mat_mul_u52(C, A, B);
                            if (ok[f])
                                TIME_BEST(tm[f], reps, nmod_mat_mul_u52(C, A, B));
                            break;
                        case FN_K52:
                            ok[f] = nmod_mat_mul_k52(C, A, B);
                            if (ok[f])
                                TIME_BEST(tm[f], reps, nmod_mat_mul_k52(C, A, B));
                            break;
                        case FN_FP50:
                            ok[f] = nmod_mat_mul_fp50(C, A, B);
                            if (ok[f])
                                TIME_BEST(tm[f], reps, nmod_mat_mul_fp50(C, A, B));
                            break;
                        case FN_BLAS:
                            ok[f] = nmod_mat_mul_blas(C, A, B);
                            if (ok[f])
                                TIME_BEST(tm[f], reps, nmod_mat_mul_blas(C, A, B));
                            break;
                        case FN_FGEMM:
                            ok[f] = nmod_mat_mul_fgemm(C, A, B);
                            if (ok[f])
                                TIME_BEST(tm[f], reps, nmod_mat_mul_fgemm(C, A, B));
                            break;
                        case FN_CLASSICAL:
                            ok[f] = 1;
                            TIME_BEST(tm[f], reps, nmod_mat_mul_classical(C, A, B));
                            break;
                        case FN_STRASSEN:
                            ok[f] = 1;
                            TIME_BEST(tm[f], reps, nmod_mat_mul_strassen(C, A, B));
                            break;
                        case FN_THREADED:
                            ok[f] = 1;
                            TIME_BEST(tm[f], reps, nmod_mat_mul_classical_threaded(C, A, B));
                            break;
                        default:
                            ok[f] = 1;
                            TIME_BEST(tm[f], reps, nmod_mat_mul(C, A, B));
                            break;
                    }

                    if (ok[f] && tm[f] > tmax)
                        skipped[f] = 1;

                    if (ok[f] && csv)
                        flint_printf("csv,%wd,%wu,%wd,%wd,%wd,%wd,%s,%.4f\n",
                                     pbits, p, nthreads[t], m, k, n,
                                     fn_names[f], tm[f] * 1e6);
                }

                /* the fastest function overall (mul excluded, being a
                   dispatch of the others), and the fastest SIMD kernel and
                   base algorithm for the crossover ratio */
                for (f = 0; f < NUM_FN; f++)
                {
                    if (!ok[f] || f == FN_MUL)
                        continue;
                    if (best_fn < 0 || tm[f] < tm[best_fn])
                        best_fn = f;
                    if (fn_is_base[f] && tm[f] < best_other)
                        best_other = tm[f];
                    if (fn_is_simd[f] && tm[f] < best_simd)
                        best_simd = tm[f];
                }

                flint_printf("%wd x %wd x %wd", m, k, n);
                shape_len = 6;
                {
                    slong tmp;
                    for (tmp = m; tmp > 0 || shape_len == 6; tmp /= 10) shape_len++;
                    for (tmp = k; tmp > 0; tmp /= 10) shape_len++;
                    for (tmp = n; tmp > 0; tmp /= 10) shape_len++;
                }
                while (shape_len++ < 18)
                    flint_printf(" ");

                for (f = 0; f < NUM_FN; f++)
                {
                    if (!do_fn[f])
                        continue;
                    if (ok[f])
                        flint_printf(" %10.1f", tm[f] * 1e6);
                    else
                        flint_printf(" %10s", "-");
                }

                if (best_fn >= 0)
                    flint_printf("   %-10s", fn_names[best_fn]);
                else
                    flint_printf("   %-10s", "-");

                if (best_simd < 1e300 && best_other < 1e300)
                {
                    flint_printf(" %.3f", best_simd / best_other);

                    last_dim_measured = FLINT_MIN(FLINT_MIN(m, k), n);
                    if (best_simd < best_other)
                    {
                        if (first_win < 0)
                            first_win = last_dim_measured;
                    }
                    else
                    {
                        if (first_win >= 0)
                            simd_ever_lost_after = 1;
                        first_win = -1;
                    }
                }

                flint_printf("\n");
                fflush(stdout);

                if (do_nmod32)
                    nmod32_data_clear(&w32);
                nmod_mat_clear(A);
                nmod_mat_clear(B);
                nmod_mat_clear(C);
            }

            if (!have_shape)
            {
                if (first_win >= 0)
                    flint_printf("  => a SIMD kernel beats blas/fgemm/classical/threaded from dim %wd on%s\n",
                                 first_win, simd_ever_lost_after ?
                                 " (after winning and losing at smaller sizes)" : "");
                else if (last_dim_measured > 0)
                    flint_printf("  => no SIMD kernel beats all of blas/fgemm/classical/threaded at dim %wd\n",
                                 last_dim_measured);
            }

            flint_printf("\n");
        }
    }

    flint_rand_clear(state);

    return 0;
}
