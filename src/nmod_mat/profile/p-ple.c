/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Timings of the PLE decomposition.

    Usage: p-ple [mode] [bits] [nmax]

      mode = square  (default) n x n uniformly random matrices: all the
                     algorithms of nmod_mat_ple, nmod_mat_lu and
                     nmod_mat_mul (n x n by n x n)
      mode = shapes  m x n uniformly random matrices, m = n/8, n/2, 2n, 8n
      mode = ranks   n x n matrices L R U, with L, U random triangular and
                     R a random rank profile matrix of rank n, n/2, n/8
                     (many upward row rotations; the rank n/2, n/8 matrices have
                     many columns without pivot)

      bits = comma-separated list of bit sizes of the prime (default:
             2,8,16,20,26,28,29,31,32,40,50,52,55,60,63,64)
      nmax = largest dimension (default 1024; 512 for shapes)

    Times in microseconds per call, including a copy of the input (not for
    mul); the classical algorithm is not timed above min(m, n) = 256, the
    delayed and Crout algorithms above 1024 (columns "-"). The ratios
    compare nmod_mat_ple with nmod_mat_lu and nmod_mat_mul.
*/

#include <stdlib.h>
#include <string.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

#define NUM_ALGOS 7
static const char * algo_name[NUM_ALGOS] =
    { "ple", "classical", "delayed", "crout", "recursive", "lu", "mul" };

/* runs the algorithm reps times on a copy of A; B, C, P, Q workspace */
static void
run_algo(int algo, slong reps, nmod_mat_t B, nmod_mat_t C,
         slong * P, slong * Q, const nmod_mat_t A)
{
    slong i;

    for (i = 0; i < reps; i++)
    {
        if (algo == 6)
        {
            nmod_mat_mul(C, A, A);
            continue;
        }

        nmod_mat_set(B, A);
        if (algo == 0)
            nmod_mat_ple(P, Q, B, 0);
        else if (algo == 1)
            nmod_mat_ple_classical(P, Q, B, 0);
        else if (algo == 2)
            nmod_mat_ple_classical_delayed(P, Q, B, 0);
        else if (algo == 3)
            nmod_mat_ple_crout(P, Q, B, 0);
        else if (algo == 4)
            nmod_mat_ple_recursive(P, Q, B, 0);
        else
            nmod_mat_lu(P, B, 0);
    }
}

/*
    Times in microseconds per call of the algorithms on A (negative: not
    timed). Each algorithm is repeated in batches of at least 10 ms, and
    the batches of the different algorithms are interleaved, five rounds,
    keeping the smallest time per algorithm: the ratios are then less
    sensitive to slow phases of the machine.
*/
static void
time_algos(double * t, const nmod_mat_t A, int with_mul)
{
    slong m = A->r, n = A->c, s = FLINT_MIN(m, n);
    slong * P = flint_malloc((m + 1) * sizeof(slong));
    slong * Q = flint_malloc((n + 1) * sizeof(slong));
    slong reps[NUM_ALGOS];
    int enabled[NUM_ALGOS];
    nmod_mat_t B, C;
    timeit_t timer;
    int algo, round;

    nmod_mat_init(B, m, n, A->mod.n);
    nmod_mat_init(C, n, n, A->mod.n);

    enabled[0] = enabled[4] = enabled[5] = 1;
    enabled[1] = (s <= 256);
    enabled[2] = enabled[3] = (s <= 1024);
    enabled[6] = (with_mul && m == n);

    for (algo = 0; algo < NUM_ALGOS; algo++)
    {
        t[algo] = -1.0;
        if (!enabled[algo])
            continue;

        for (reps[algo] = 1; ; reps[algo] *= 2)
        {
            timeit_start_us(timer);
            run_algo(algo, reps[algo], B, C, P, Q, A);
            timeit_stop_us(timer);
            if (timer->wall >= 10000)
                break;
        }
    }

    for (round = 0; round < 5; round++)
    {
        for (algo = 0; algo < NUM_ALGOS; algo++)
        {
            double tt;

            if (!enabled[algo])
                continue;

            timeit_start_us(timer);
            run_algo(algo, reps[algo], B, C, P, Q, A);
            timeit_stop_us(timer);

            tt = (double) timer->wall / reps[algo];
            if (t[algo] < 0 || tt < t[algo])
                t[algo] = tt;
        }
    }

    nmod_mat_clear(B);
    nmod_mat_clear(C);
    flint_free(P);
    flint_free(Q);
}

static void
print_header(const char * first)
{
    int algo;

    flint_printf("%s", first);
    for (algo = 0; algo < NUM_ALGOS; algo++)
        flint_printf(" %11s", algo_name[algo]);
    flint_printf("   ple/lu ple/mul\n");
}

static void
print_times(const double * t)
{
    int algo;

    for (algo = 0; algo < NUM_ALGOS; algo++)
    {
        if (t[algo] < 0)
            flint_printf(" %11s", "-");
        else
            flint_printf(" %11.2f", t[algo]);
    }
    flint_printf("   %6.3f", t[0] / t[5]);
    if (t[6] > 0)
        flint_printf("  %6.3f", t[0] / t[6]);
    flint_printf("\n");
}

/* A = L R U, L (m x m) and U (n x n) random unit triangular, R with
   random nonzero entries at r random positions in distinct rows and
   columns */
static void
randrank_profile_dense(nmod_mat_t A, slong r, flint_rand_t state)
{
    slong m = A->r, n = A->c, k, i;
    slong * mrp_r = flint_malloc((m + 1) * sizeof(slong));
    slong * mrp_c = flint_malloc((n + 1) * sizeof(slong));
    nmod_mat_t L, U;

    for (k = 0; k < m; k++)
        mrp_r[k] = k;
    for (k = 0; k < n; k++)
        mrp_c[k] = k;
    for (k = 0; k < r; k++)
    {
        FLINT_SWAP(slong, mrp_r[k], mrp_r[k + n_randint(state, m - k)]);
        FLINT_SWAP(slong, mrp_c[k], mrp_c[k + n_randint(state, n - k)]);
    }

    nmod_mat_init(L, m, r, A->mod.n);
    nmod_mat_init(U, r, n, A->mod.n);
    nmod_mat_rand(L, state);
    nmod_mat_rand(U, state);
    for (k = 0; k < r; k++)
    {
        for (i = 0; i < mrp_r[k]; i++)
            nmod_mat_entry(L, i, k) = 0;
        nmod_mat_entry(L, mrp_r[k], k) = 1;
        for (i = 0; i < mrp_c[k]; i++)
            nmod_mat_entry(U, k, i) = 0;
        nmod_mat_entry(U, k, mrp_c[k]) = 1 + n_randint(state, A->mod.n - 1);
    }
    nmod_mat_mul(A, L, U);

    nmod_mat_clear(L);
    nmod_mat_clear(U);
    flint_free(mrp_r);
    flint_free(mrp_c);
}

static const slong default_bits[] =
    { 2, 8, 16, 20, 26, 28, 29, 31, 32, 40, 50, 52, 55, 60, 63, 64, 0 };

static const slong square_sizes[] = { 4, 8, 12, 16, 24, 32, 48, 64, 96,
    128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 0 };

/* a prime of the given bit size: 2 for 1 or 2 bits, the largest prime
   below 2^FLINT_BITS for FLINT_BITS bits, otherwise the smallest prime
   above 2^(bits - 1) */
static ulong
prime_of_bits(slong bits)
{
    if (bits <= 2)
        return 2;
    if (bits >= FLINT_BITS)
#if FLINT_BITS == 64
        return UWORD(18446744073709551557);
#else
        return UWORD(4294967291);
#endif
    return n_nextprime(UWORD(1) << (bits - 1), 0);
}

int main(int argc, char ** argv)
{
    const char * mode = (argc > 1) ? argv[1] : "square";
    slong bits_list[65], nb = 0, nmax, ib, n;
    flint_rand_t state;
    double t[NUM_ALGOS];

    if (argc > 2)
    {
        char * s = argv[2], * e;
        while (*s != '\0' && nb < 64)
        {
            slong b = strtol(s, &e, 10);
            if (e == s || b < 1 || b > FLINT_BITS)
            {
                flint_printf("invalid list of bit sizes: %s\n", argv[2]);
                return 1;
            }
            bits_list[nb++] = b;
            s = (*e == ',') ? e + 1 : e;
        }
    }
    else
    {
        for (nb = 0; default_bits[nb] != 0; nb++)
            bits_list[nb] = default_bits[nb];
    }

    nmax = (argc > 3) ? atol(argv[3]) : (strcmp(mode, "shapes") ? 1024 : 512);

    flint_rand_init(state);
    flint_set_num_threads(1);

    if (strcmp(mode, "square") == 0)
    {
        flint_printf("n x n uniformly random, times in us\n");
        print_header("bits      n");
        for (ib = 0; ib < nb; ib++)
        {
            ulong p = prime_of_bits(bits_list[ib]);

            slong is;

            for (is = 0; (n = square_sizes[is]) != 0 && n <= nmax; is++)
            {
                nmod_mat_t A;

                nmod_mat_init(A, n, n, p);
                nmod_mat_rand(A, state);
                time_algos(t, A, 1);
                flint_printf("%4wd %6wd", bits_list[ib], n);
                print_times(t);
                nmod_mat_clear(A);
            }
        }
    }
    else if (strcmp(mode, "shapes") == 0)
    {
        flint_printf("m x n uniformly random, times in us\n");
        print_header("bits      m      n");
        for (ib = 0; ib < nb; ib++)
        {
            ulong p = prime_of_bits(bits_list[ib]);

            for (n = 32; n <= nmax; n *= 2)
            {
                slong shape;
                slong ms[4] = { n / 8, n / 2, 2 * n, 8 * n };

                for (shape = 0; shape < 4; shape++)
                {
                    nmod_mat_t A;

                    nmod_mat_init(A, ms[shape], n, p);
                    nmod_mat_rand(A, state);
                    time_algos(t, A, 0);
                    flint_printf("%4wd %6wd %6wd", bits_list[ib], ms[shape], n);
                    print_times(t);
                    nmod_mat_clear(A);
                }
            }
        }
    }
    else if (strcmp(mode, "ranks") == 0)
    {
        flint_printf("n x n, L R U with a random rank profile matrix R of rank r, times in us\n");
        print_header("bits      n      r");
        for (ib = 0; ib < nb; ib++)
        {
            ulong p = prime_of_bits(bits_list[ib]);

            for (n = 16; n <= nmax; n *= 2)
            {
                slong ir;
                slong rs[3] = { n, n / 2, n / 8 };

                for (ir = 0; ir < 3; ir++)
                {
                    nmod_mat_t A;

                    nmod_mat_init(A, n, n, p);
                    randrank_profile_dense(A, rs[ir], state);
                    time_algos(t, A, 1);
                    flint_printf("%4wd %6wd %6wd", bits_list[ib], n, rs[ir]);
                    print_times(t);
                    nmod_mat_clear(A);
                }
            }
        }
    }
    else
    {
        flint_printf("unknown mode %s (square, shapes, ranks)\n", mode);
        flint_rand_clear(state);
        return 1;
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
