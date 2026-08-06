/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Tuning program for the quadratic sieve.

    Usage: tune-qsieve BITS [-num N] [-reps REPS] [-delta DELTA] [-threads THREADS]

      BITS     bit size of the semiprimes to tune for (required)
      N        number of random semiprimes to factor per evaluation (default 100)
      REPS     times each evaluation is repeated, taking the fastest (default 1)
      DELTA    search parameter (default 0.05)
      THREADS  number of threads to profile (the wall time is counted)

    The quadratic sieve has five tuning parameters, currently hardcoded per
    bit-size in the qsieve_tune table in qsieve.h:

       ks_primes     number of Knuth-Schroeppel primes
       fb_primes     number of factor base primes (including k and 2)
       small_primes  number of small primes not sieved with
       sieve_size    size of the sieve interval
       sieve_bits    sieve threshold

    The objective (total time to factor all N semiprimes) is non-convex, noisy,
    and can fail to converge for extreme parameter values.  

    We use a robust coordinate-descent local search: starting from the table
    defaults in qsieve.h, each round we perturb every parameter in turn, up
    and down by the relative factor DELTA (or at least by an absolute step size
    of 1), evaluate each single-parameter change, and adopt the best one if it
    beats the (freshly re-measured) baseline.
    We stop when no perturbation improves.

    The search is likely to find a local minimum rather than a global one,
    so the initial table parameters should be reasonable to begin with.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_factor.h"
#include "ulong_extras.h"
#include "qsieve.h"

#define NPARAMS 5

static const char * param_name[NPARAMS] =
{
    "ks_primes", "fb_primes", "small_primes", "sieve_size", "sieve_bits"
};

static double
wall_clock(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
}

/* generate a balanced semiprime of exactly `bits` bits */
static void
make_semiprime(fmpz_t n, flint_rand_t state, slong bits)
{
    slong hb = bits / 2;
    fmpz_t p, q;
    fmpz_init(p);
    fmpz_init(q);

    do
    {
        fmpz_randbits(p, state, hb);
        fmpz_abs(p, p);
        fmpz_nextprime(p, p, 1);
        fmpz_randbits(q, state, bits - hb);
        fmpz_abs(q, q);
        fmpz_nextprime(q, q, 1);
        fmpz_mul(n, p, q);
    } while (fmpz_bits(n) != (flint_bitcnt_t) bits || fmpz_equal(p, q));

    fmpz_clear(p);
    fmpz_clear(q);
}

/* clamp parameter `idx` to a sane range given the other parameters in p */
static slong
clamp_param(slong idx, slong val, const slong p[NPARAMS])
{
    slong lo = 1, hi = 1;

    switch (idx)
    {
        case 0: lo = 10;            hi = 100000;       break; /* ks_primes    */
        case 1:                                               /* fb_primes    */
                lo = p[2] + 1;
                if (lo < 10)
                    lo = 10;
                hi = 2000000;
                break;
        case 2: lo = 3;             hi = p[1] - 1;     break; /* small_primes */
        case 3: lo = 1000;          hi = 64 * 65536;   break; /* sieve_size   */
        case 4: lo = 16;            hi = 250;          break; /* sieve_bits   */
    }

    if (val < lo) val = lo;
    if (val > hi) val = hi;
    return val;
}

/* total time (seconds) to factor all N semiprimes with parameters p,
   taking the fastest of `reps` repetitions; verifies correctness */
static double
eval_params(const fmpz * semis, slong N, const slong p[NPARAMS], int reps)
{
    double best = -1.0;
    int r;
    slong i;

    for (r = 0; r < reps; r++)
    {
        double t0 = wall_clock();

        for (i = 0; i < N; i++)
        {
            fmpz_factor_t fac;
            fmpz_factor_init(fac);
            qsieve_factor_with_tune(fac, semis + i,
                (ulong) p[0], p[1], p[2], p[3], (ulong) p[4]);

            /* sanity check: product of factors must equal the input */
            {
                fmpz_t chk, pw;
                slong j;
                fmpz_init_set_ui(chk, 1);
                fmpz_init(pw);
                for (j = 0; j < fac->num; j++)
                {
                    fmpz_pow_ui(pw, fac->p + j, fac->exp[j]);
                    fmpz_mul(chk, chk, pw);
                }
                if (!fmpz_equal(chk, semis + i))
                    flint_printf("    WARNING: incorrect factorization!\n");
                fmpz_clear(chk);
                fmpz_clear(pw);
            }

            fmpz_factor_clear(fac);
        }

        {
            double t = wall_clock() - t0;
            if (best < 0.0 || t < best)
                best = t;
        }
    }

    return best;
}

static void
print_params(const slong p[NPARAMS])
{
    flint_printf("[ks=%wd fb=%wd small=%wd sieve_size=%wd sieve_bits=%wd]",
                 p[0], p[1], p[2], p[3], p[4]);
}

/* Just to make the table prettier */
slong round_param(slong n)
{
    if (n < 100)
    {

    }
    else if (n < 1000)
    {
        n = (slong) (n * 0.1 + 0.5) * 10;
    }
    else if (n < 10000)
    {
        n = (slong) (n * 0.01 + 0.5) * 100;
    }
    else
    {
        n = (slong) (n * 0.001 + 0.5) * 1000;
    }

    return n;
}

int
main(int argc, char ** argv)
{
    slong bits, N, i, num_threads;
    flint_rand_t state;
    fmpz * semis;
    slong cur[NPARAMS];
    double base_time, best_overall;
    double DELTA;
    slong reps, round;

    if (argc < 2)
    {
        flint_printf("usage: %s BITS [-num N] [-reps REPS] [-delta DELTA] [-threads THREADS]\n", argv[0]);
        flint_printf("  tunes the quadratic sieve parameters for BITS-bit semiprimes\n");
        return 1;
    }

    num_threads = 1;
    N = 100;
    reps = 1;
    DELTA = 0.05;

    bits = atol(argv[1]);

    for (i = 2; i < argc; i++)
    {
        if (!strcmp(argv[i], "-threads"))
        {
            num_threads = atol(argv[i+1]);
            flint_set_num_threads(num_threads);
            i++;
        }
        else if (!strcmp(argv[i], "-num"))
        {
            N = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-reps"))
        {
            reps = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-delta"))
        {
            reps = atof(argv[i+1]);
            i++;
        }
    }

    if (bits < 10)
    {
        flint_printf("BITS should be at least 10 (quadratic sieve is for larger n)\n");
        return 1;
    }
    if (N < 1) N = 1;
    if (reps < 1) reps = 1;

    flint_rand_init(state);

    /* ---- generate the fixed set of test semiprimes ---- */
    flint_printf("Tuning quadratic sieve for %wd-bit semiprimes and %wd threads\n", bits, num_threads);
    flint_printf("Generating %wd random balanced semiprime(s), reps=%wd per evaluation\n\n", N, reps);

    semis = _fmpz_vec_init(N);
    for (i = 0; i < N; i++)
        make_semiprime(semis + i, state, bits);

    /* ---- look up the default tuning parameters for this bit-size ---- */
    {
        slong idx;
        for (idx = 1; idx < (slong) QS_TUNE_SIZE; idx++)
            if (qsieve_tune[idx][0] > (ulong) bits)
                break;
        idx--;

        cur[0] = (slong) qsieve_tune[idx][1]; /* ks_primes    */
        cur[1] = (slong) qsieve_tune[idx][2]; /* fb_primes    */
        cur[2] = (slong) qsieve_tune[idx][3]; /* small_primes */
        cur[3] = (slong) qsieve_tune[idx][4]; /* sieve_size   */
        cur[4] = (slong) qsieve_tune[idx][5]; /* sieve_bits   */

        cur[0] = round_param(cur[0]);
        cur[1] = round_param(cur[1]);
        cur[2] = round_param(cur[2]);
        cur[3] = round_param(cur[3]);
        cur[4] = round_param(cur[4]);
    }

    base_time = eval_params(semis, N, cur, reps);
    best_overall = base_time;
    flint_printf("Default parameters: ");
    print_params(cur);
    flint_printf("\n   baseline time: %.4f s (%.2f ms/semiprime)\n\n",
                 base_time, 1000.0 * base_time / (double) N);

    /* ---- coordinate-descent local search ---- */
    round = 0;
    while (1)
    {
        slong best_params[NPARAMS];
        double round_best;
        int improved = 0;
        slong k;
        int dir;

        round++;
        flint_printf("=== Round %wd ===\n", round);

        /* re-measure the current parameters as this round's baseline (the
           objective is noisy, so we always compare against a fresh reading) */
        base_time = eval_params(semis, N, cur, reps);
        round_best = base_time;
        for (k = 0; k < NPARAMS; k++)
            best_params[k] = cur[k];
        flint_printf("  baseline ");
        print_params(cur);
        flint_printf(" -> %.4f s\n", base_time);

        for (k = 0; k < NPARAMS; k++)
        {
            for (dir = 0; dir < 2; dir++)
            {
                slong cand[NPARAMS];
                slong j, v, delta, newv;
                double t;

                for (j = 0; j < NPARAMS; j++)
                    cand[j] = cur[j];

                v = cur[k];

                delta = (slong) (DELTA * (double) v + 0.5);
                if (delta < 1)
                    delta = 1;
                newv = (dir == 0) ? v + delta : v - delta;

                newv = round_param(newv);

                cand[k] = clamp_param(k, newv, cand);

                if (cand[k] == cur[k])
                {
                    flint_printf("    %-12s %s: no change after clamping (skip)\n",
                                 param_name[k], dir == 0 ? "up  " : "down");
                    continue;
                }

                t = eval_params(semis, N, cand, reps);

                flint_printf("    %-12s %s: %wd -> %wd : %.4f s  (%+.1f%%)%s\n",
                             param_name[k], dir == 0 ? "up  " : "down",
                             cur[k], cand[k], t,
                             100.0 * (t - base_time) / base_time,
                             (t < round_best) ? "  * best so far" : "");

                if (t < round_best)
                {
                    round_best = t;
                    for (j = 0; j < NPARAMS; j++)
                        best_params[j] = cand[j];
                    improved = 1;
                }
            }
        }

        if (improved)
        {
            for (k = 0; k < NPARAMS; k++)
                cur[k] = best_params[k];
            best_overall = round_best;
            flint_printf("  -> adopted ");
            print_params(cur);
            flint_printf(" : %.4f s (%.1f%% faster than round baseline)\n\n",
                         round_best, 100.0 * (base_time - round_best) / base_time);
        }
        else
        {
            flint_printf("  -> no improvement; stopping.\n\n");
            break;
        }
    }

    flint_printf("================ RESULT ================\n");
    flint_printf("Best parameters for %wd bits: ", bits);
    print_params(cur);
    flint_printf("\n");
    flint_printf("Tuning-table row format:     {%wd, %wd, %wd, %wd, %wd, %wd},  /* %wd digits, %.2f ms/semiprime */ \n",
                 bits, cur[0], cur[1], cur[2], cur[3], cur[4],
                (slong) (bits * 0.30103), 1000.0 * best_overall / (double) N);
    flint_printf("Final time: %.4f s (%.2f ms/semiprime)\n",
                 best_overall, 1000.0 * best_overall / (double) N);

    _fmpz_vec_clear(semis, N);
    flint_rand_clear(state);
    flint_cleanup();
    return 0;
}
