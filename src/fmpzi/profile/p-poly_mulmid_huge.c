/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Time and memory of the transformed algorithms for large Gaussian
    integer polynomial products.

    Usage: p-poly_mulmid_huge [options]

      -s SHAPE   mul  n x n -> [0, 2n - 1)            (default)
                 sqr  squaring
                 mid  (2n - 1) x n -> [n - 1, 2n - 1)
      -b B       smallest part size in bits (default 16384); sizes
                 B, 4B, 16B, ... up to 2^20
      -w W       largest n * bits (default 3e8)
      -W W       smallest n * bits (default 1e6)
      -m MODE    both (default), time or mem
      -t T       number of threads (default 1)

    For each size (bits per real and imaginary part, length n) it prints
    the extra peak memory in MB, then (in a second pass) the time in ms,
    of

      KS      _fmpzi_poly_mulmid_KS
      pair    _fmpzi_poly_mulmid_SS_pair
      single  _fmpzi_poly_mulmid_SS_single
      disp    the dispatcher _fmpzi_poly_mulmid

    The memory is measured on Linux only, each call in a child process
    forked from a parent that has not computed anything yet (one that
    had would pass on warm buffers): the peak resident size
    during the call minus the resident size before it, the inputs and
    outputs (preallocated at their final size) excluded.
*/

#include <stdlib.h>
#include <string.h>
#include "longlong.h"
#include "fmpz.h"
#include "fmpzi.h"
#include "thread_support.h"
#include "profiler.h"

#if defined(__linux__)
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/resource.h>
#define HAVE_MEM 1
#else
#define HAVE_MEM 0
#endif

typedef void (*mulmid_fn)(fmpzi_struct *, const fmpzi_struct *, slong,
                          const fmpzi_struct *, slong, slong, slong);

#define NALG 4
static const mulmid_fn algs[NALG] = { _fmpzi_poly_mulmid_KS,
    _fmpzi_poly_mulmid_SS_pair, _fmpzi_poly_mulmid_SS_single, _fmpzi_poly_mulmid };
static const char * names[NALG] = { "KS", "pair", "single", "disp" };

typedef struct
{
    fmpzi_struct * a, * b, * r;
    slong len1, len2, nlo, nhi;
}
prob_struct;

static void
prob_init(prob_struct * P, int shape, slong n, slong bits, flint_rand_t state)
{
    slong i, w;

    P->len2 = n;
    P->len1 = (shape == 2) ? 2 * n - 1 : n;
    P->nlo = (shape == 2) ? n - 1 : 0;
    P->nhi = (shape == 2) ? 2 * n - 1 : P->len1 + P->len2 - 1;
    w = P->nhi - P->nlo;

    P->a = flint_malloc(P->len1 * sizeof(fmpzi_struct));
    for (i = 0; i < P->len1; i++)
    {
        fmpzi_init(P->a + i);
        fmpz_randbits(fmpzi_realref(P->a + i), state, bits);
        fmpz_randbits(fmpzi_imagref(P->a + i), state, bits);
    }

    if (shape == 1)
        P->b = P->a;
    else
    {
        P->b = flint_malloc(P->len2 * sizeof(fmpzi_struct));
        for (i = 0; i < P->len2; i++)
        {
            fmpzi_init(P->b + i);
            fmpz_randbits(fmpzi_realref(P->b + i), state, bits);
            fmpz_randbits(fmpzi_imagref(P->b + i), state, bits);
        }
    }

    /* outputs at their final size, so that the measurement does not
       include them */
    P->r = flint_malloc(w * sizeof(fmpzi_struct));
    for (i = 0; i < w; i++)
    {
        fmpzi_init(P->r + i);
        fmpz_one(fmpzi_realref(P->r + i));
        fmpz_mul_2exp(fmpzi_realref(P->r + i), fmpzi_realref(P->r + i),
                      2 * bits + FLINT_BIT_COUNT(n) + 2);
        fmpz_set(fmpzi_imagref(P->r + i), fmpzi_realref(P->r + i));
    }
}

static void
prob_clear(prob_struct * P, int shape)
{
    slong i;

    for (i = 0; i < P->len1; i++)
        fmpzi_clear(P->a + i);
    flint_free(P->a);
    if (shape != 1)
    {
        for (i = 0; i < P->len2; i++)
            fmpzi_clear(P->b + i);
        flint_free(P->b);
    }
    for (i = 0; i < P->nhi - P->nlo; i++)
        fmpzi_clear(P->r + i);
    flint_free(P->r);
}

static void
prob_run(prob_struct * P, int alg)
{
    algs[alg](P->r, P->a, P->len1, P->b, P->len2, P->nlo, P->nhi);
}

/* best of up to five batches of at least 20 ms; a single call when
   that takes more than 2 s */
static double
time_ms(prob_struct * P, int alg)
{
    double best = 1e300;
    slong reps, i, k;
    timeit_t t;

    for (k = 0; k < 5; k++)
    {
        reps = 1;
        for (;;)
        {
            timeit_start_us(t);
            for (i = 0; i < reps; i++)
                prob_run(P, alg);
            timeit_stop_us(t);
            if (t->wall >= 20000)
                break;
            reps *= (t->wall < 1000) ? 10 : 2;
        }
        best = FLINT_MIN(best, (double) t->wall / reps);
        if (t->wall > 2000000)
            break;
    }

    return best / 1000.0;
}

#if HAVE_MEM
static double
rss_mb(void)
{
    long pages = 0, res = 0;
    FILE * f = fopen("/proc/self/statm", "r");

    if (f == NULL)
        return 0.0;
    if (fscanf(f, "%ld %ld", &pages, &res) != 2)
        res = 0;
    fclose(f);
    return (double) res * (double) sysconf(_SC_PAGESIZE) / 1048576.0;
}

static double
mem_mb(int alg, int shape, slong n, slong bits, slong threads)
{
    int fd[2];
    pid_t pid;
    double v = -1.0;

    if (pipe(fd) != 0)
        return -1.0;

    pid = fork();
    if (pid == 0)
    {
        prob_struct P;
        struct rusage ru;
        double base;
        flint_rand_t state;

        flint_rand_init(state);
        flint_set_num_threads(threads);
        prob_init(&P, shape, n, bits, state);
        base = rss_mb();
        prob_run(&P, alg);
        getrusage(RUSAGE_SELF, &ru);
        v = (double) ru.ru_maxrss / 1024.0 - base;
        if (write(fd[1], &v, sizeof(v)) != sizeof(v))
            _exit(1);
        _exit(0);
    }

    if (pid < 0 || read(fd[0], &v, sizeof(v)) != sizeof(v))
        v = -1.0;
    if (pid > 0)
        waitpid(pid, NULL, 0);
    close(fd[0]);
    close(fd[1]);
    return v;
}
#endif

int main(int argc, char * argv[])
{
    static const char * shapes[] = { "mul n x n", "sqr n", "mid (2n - 1) x n" };
    slong minbits = 16384, threads = 1, n, bits;
    double maxwork = 3e8, minwork = 1e6;
    int shape = 0, do_time = 1, do_mem = HAVE_MEM, i, alg, pass;
    flint_rand_t state;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-s") && i + 1 < argc)
        {
            i++;
            shape = !strcmp(argv[i], "sqr") ? 1 : !strcmp(argv[i], "mid") ? 2 : 0;
        }
        else if (!strcmp(argv[i], "-b") && i + 1 < argc)
        {
            minbits = atol(argv[++i]);
            minbits = FLINT_MAX(minbits, 1);
        }
        else if (!strcmp(argv[i], "-w") && i + 1 < argc)
            maxwork = atof(argv[++i]);
        else if (!strcmp(argv[i], "-W") && i + 1 < argc)
            minwork = atof(argv[++i]);
        else if (!strcmp(argv[i], "-t") && i + 1 < argc)
        {
            threads = atol(argv[++i]);
            threads = FLINT_MAX(threads, 1);
        }
        else if (!strcmp(argv[i], "-m") && i + 1 < argc)
        {
            i++;
            do_time = strcmp(argv[i], "mem") != 0;
            do_mem = HAVE_MEM && strcmp(argv[i], "time") != 0;
        }
        else
        {
            flint_printf("usage: %s [-s mul|sqr|mid] [-b minbits] [-w maxwork]"
                         " [-W minwork] [-m both|time|mem] [-t threads]\n", argv[0]);
            return 1;
        }
    }

    flint_rand_init(state);

    /* two passes over the sizes, memory first: the measuring children
       must be forked from a process that has not computed anything */
    for (pass = 0; pass < 2; pass++)
    {
        if ((pass == 0 && !do_mem) || (pass == 1 && !do_time))
            continue;

        flint_printf("\nshape %s, threads %wd: %s\n", shapes[shape], threads,
                     pass == 0 ? "extra peak memory (MB)" : "time (ms)");
        flint_printf("%8s %7s ", "bits", "n");
        for (alg = 0; alg < NALG; alg++)
            flint_printf(" %9s", names[alg]);
        flint_printf("\n");

        for (bits = minbits; bits <= (WORD(1) << 20); bits *= 4)
        {
            for (n = 16; n <= (WORD(1) << 20); n *= 4)
            {
                double v[NALG];
                prob_struct P;

                if ((double) n * bits > maxwork || (double) n * bits < minwork)
                    continue;

                if (pass == 1)
                {
                    flint_set_num_threads(threads);
                    prob_init(&P, shape, n, bits, state);
                }

                for (alg = 0; alg < NALG; alg++)
                {
#if HAVE_MEM
                    if (pass == 0)
                        v[alg] = mem_mb(alg, shape, n, bits, threads);
                    else
#endif
                        v[alg] = time_ms(&P, alg);
                }

                if (pass == 1)
                    prob_clear(&P, shape);

                flint_printf("%8wd %7wd ", bits, n);
                for (alg = 0; alg < NALG; alg++)
                {
                    if (v[alg] < 0)
                        flint_printf(" %9s", "-");
                    else
                        flint_printf(pass == 0 ? " %9.1f" : " %9.2f", v[alg]);
                }
                flint_printf("\n");
                fflush(stdout);
            }
        }
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
