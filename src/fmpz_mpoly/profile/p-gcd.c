/*
    Copyright 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* usage
likwid-setFrequencies -c S0 -g powersave
likwid-setFrequencies -g performance
make profile MOD=fmpz_mpoly && ./build/fmpz_mpoly/profile/p-gcd 4 sparse 7 5 4 8

p-gcd nthreads sparse m1 n1 m2 n2:
    run the sparse benchmark on nthreads with powers (m1, n1) and (m2, n2)
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "fmpz_mpoly.h"

#define CALCULATE_MACHINE_EFFICIENCY 0

int * cpu_affinities;

#if CALCULATE_MACHINE_EFFICIENCY

typedef struct _worker_arg_struct
{
    fmpz_mpoly_t G;
    const fmpz_mpoly_struct * A, * B;
    const fmpz_mpoly_ctx_struct * ctx;
} worker_arg_struct;

typedef worker_arg_struct worker_arg_t[1];


static void worker_gcd(void * varg)
{
    worker_arg_struct * W = (worker_arg_struct *) varg;
    fmpz_mpoly_gcd_threaded(W->G, W->A, W->B, W->ctx, 1);
}

#endif

void profile_gcd(
    const fmpz_mpoly_t realG,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    slong max_threads)
{
    fmpz_mpoly_t G;
    timeit_t timer;
    slong num_threads;
    slong serial_time;

    flint_set_num_threads(1);
    flint_set_thread_affinity(cpu_affinities, 1);
    fmpz_mpoly_init(G, ctx);
    timeit_start(timer);
    fmpz_mpoly_gcd(G, A, B, ctx);
    timeit_stop(timer);
    serial_time = FLINT_MAX(WORD(1), timer->wall);
    flint_printf("serial time: %wd\n", serial_time);
    if (!fmpz_mpoly_equal(G, realG, ctx))
    {
        printf("gcd wrong\n");
        flint_abort();
    }

    for (num_threads = 2; num_threads <= max_threads; num_threads++)
    {
        slong parallel_time;
        double parallel_efficiency;
#if CALCULATE_MACHINE_EFFICIENCY
        thread_pool_handle * handles;
        worker_arg_struct * worker_args;
        slong i;
        double machine_efficiency;
        slong num_workers;
#endif

        flint_set_num_threads(num_threads);
        flint_set_thread_affinity(cpu_affinities, num_threads);

#if CALCULATE_MACHINE_EFFICIENCY
        handles = (thread_pool_handle *) flint_malloc((num_threads - 1)*sizeof(thread_pool_handle));
        num_workers = thread_pool_request(global_thread_pool, handles, num_threads - 1);
        worker_args = (worker_arg_struct *) flint_malloc((num_workers + 1)*sizeof(worker_arg_t));

        timeit_start(timer);
        for (i = 0; i <= num_workers; i++)
        {
            fmpz_mpoly_init((worker_args + i)->G, ctx);
            (worker_args + i)->A = A;
            (worker_args + i)->B = B;
            (worker_args + i)->ctx = ctx;
            if (i < num_workers)
            {
                thread_pool_wake(global_thread_pool, handles[i], 0, worker_gcd, worker_args + i);
            }
            else
            {
                worker_gcd(worker_args + i);
            }
        }
        for (i = 0; i < num_workers; i++)
        {
            thread_pool_wait(global_thread_pool, handles[i]);
        }
        timeit_stop(timer);
        parallel_time = FLINT_MAX(WORD(1), timer->wall);

        for (i = 0; i <= num_workers; i++)
        {
            if (!fmpz_mpoly_equal((worker_args + i)->G, realG, ctx))
            {
                printf("gcd wrong\n");
                flint_abort();
            }
            fmpz_mpoly_clear((worker_args + i)->G, ctx);

            if (i < num_workers)
            {
                thread_pool_give_back(global_thread_pool, handles[i]);
            }
        }
        flint_free(worker_args);
        flint_free(handles);

        machine_efficiency = (double)(serial_time)/(double)(parallel_time);
#endif

        fmpz_mpoly_clear(G, ctx);
        fmpz_mpoly_init(G, ctx);
        timeit_start(timer);
        fmpz_mpoly_gcd(G, A, B, ctx);
        timeit_stop(timer);
        parallel_time = FLINT_MAX(WORD(1), timer->wall);
        if (!fmpz_mpoly_equal(G, realG, ctx))
        {
            printf("gcd wrong\n");
            flint_abort();
        }

        parallel_efficiency = (double)(serial_time)/(double)(parallel_time)/(double)(num_threads);

#if CALCULATE_MACHINE_EFFICIENCY
        flint_printf("parallel %wd time: %wd, efficiency %f (machine %f)\n", num_threads, parallel_time, parallel_efficiency, machine_efficiency);
#else
        flint_printf("parallel %wd time: %wd, efficiency %f\n", num_threads, parallel_time, parallel_efficiency);
#endif
    }

    fmpz_mpoly_clear(G, ctx);
}


void profile_power(const char * astr, const char * bstr, slong nvars,
  const char * name, slong m1, slong n1, slong m2, slong n2, slong max_threads)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t a, b, t, A, B, G;
    const char * vars[] = {"x", "y", "z", "t", "u", "v" , "w"};

    FLINT_ASSERT(nvars <= 7);

    fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);
    fmpz_mpoly_init(a, ctx);
    fmpz_mpoly_init(b, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(A, ctx);
    fmpz_mpoly_init(B, ctx);
    fmpz_mpoly_init(G, ctx);

    fmpz_mpoly_set_str_pretty(a, astr, vars, ctx);
    fmpz_mpoly_set_str_pretty(b, bstr, vars, ctx);
    fmpz_mpoly_pow_ui(A, a, m1, ctx);
    fmpz_mpoly_pow_ui(t, b, n1, ctx);
    fmpz_mpoly_mul(A, A, t, ctx);
    fmpz_mpoly_pow_ui(B, a, m2, ctx);
    fmpz_mpoly_pow_ui(t, b, n2, ctx);
    fmpz_mpoly_mul(B, B, t, ctx);
    fmpz_mpoly_pow_ui(G, a, FLINT_MIN(m1, m2), ctx);
    fmpz_mpoly_pow_ui(t, b, FLINT_MIN(n1, n2), ctx);
    fmpz_mpoly_mul(G, G, t, ctx);
    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

    flint_printf("starting %s gcd (%wu, %wd), (%wd, %wd):\n", name, m1, n1, m2, n2);
    profile_gcd(G, A, B, ctx, max_threads);

    fmpz_mpoly_clear(G, ctx);
    fmpz_mpoly_clear(B, ctx);
    fmpz_mpoly_clear(A, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(b, ctx);
    fmpz_mpoly_clear(a, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main(int argc, char *argv[])
{
    slong i, m1, n1, m2, n2, max_threads;
    const slong thread_limit = 64;
    const char * name;

    cpu_affinities = flint_malloc(thread_limit*sizeof(int));
    for (i = 0; i < thread_limit; i++)
        cpu_affinities[i] = i;

    if (argc == 7)
    {
        max_threads = atoi(argv[1]);
        max_threads = FLINT_MIN(max_threads, thread_limit);
        max_threads = FLINT_MAX(max_threads, WORD(1));
        name = argv[2];
        m1 = atoi(argv[3]);
        n1 = atoi(argv[4]);
        m2 = atoi(argv[5]);
        n2 = atoi(argv[6]);
    }
    else
    {
        printf("  usage: p-gcd nthreads sparse m1 n1 m2 n2\n");
        printf("running: p-gcd 4 sparse 7 5 4 8 \n");
        max_threads = 4;
        name = "sparse";
        m1 = 7;
        n1 = 5;
        m2 = 4;
        n2 = 8;
    }

    flint_printf("setting up fmpz_mpoly %s gcd ... ", name);

    if (strcmp(name, "dense2") == 0)
    {
        profile_power("1 + x + 2*y + x*y",
                      "1 + 2*x + y + x*y", 2, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "semisparse2") == 0)
    {
        profile_power("1 + x + 2*y",
                      "1 + x + y", 2, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "dense3") == 0)
    {
        profile_power("1 + x + y + z + x*y + x*z",
                      "1 + x + y + z + y*z + y*z", 3, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "semisparse3") == 0)
    {
        profile_power("1 + 2*x + 3*y - z",
                      "1 + x + y + z", 3, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "sparse3") == 0)
    {
        profile_power("1 + x^2 + y^20 + z^33",
                      "1 + x^21 + y^11 + z^7", 3, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "semisparse4") == 0)
    {
        profile_power("1 + 2*x + 3*y - z - 2*t",
                      "1 + x + y + z + t", 4, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "boldsparse4") == 0)
    {
        profile_power("1 + x + y^2 + z^3 + t^4",
                      "1 + x^4 + y^3 + z^2 + t", 4, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "semisparse5") == 0)
    {
        profile_power("1 + 2*x + 3*y - z - 2*t - 3*u",
                      "1 + x + y + z + t + u", 5, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "boldsparse5") == 0)
    {
        profile_power("1 + x^1 + y^2 + z^3 + t^4 + u^5",
                      "1 + x^5 + y^4 + z^3 + t^2 + u^1", 5, name, m1, n1, m2, n2, max_threads);
    }
    else if (strcmp(name, "heavysparse5") == 0)
    {
        profile_power("1 + x^1 + y^2 + z^6 + t^4 + u^9",
                      "1 + x^9 + y^4 + z^3 + t^7 + u^1", 5, name, m1, n1, m2, n2, max_threads);
    }
    else
    {
        profile_power("1 + x^1 + y^5 + z^4 + t^40 + u^50",
                      "1 + x^9 + y^2 + z^11 + t^7 + u^27", 5, "sparse", m1, n1, m2, n2, max_threads);
    }

    flint_free(cpu_affinities);

    flint_cleanup_master();
    return 0;
}

