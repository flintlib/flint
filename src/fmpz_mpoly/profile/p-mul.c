/*
    Copyright 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* usage:
likwid-setFrequencies -g performance
make profile MOD=fmpz_mpoly && ./build/fmpz_mpoly/profile/p-mul 4 sparse 12 12

p-mul nthreads sparse m n:
    run the sparse benchmark on nthreads with powers (m, n)
    mul((1+x+y+2*z^2+3*t^3+5*u^5)^m, (1+u+t+2*z^2+3*y^3+5*x^5)^n)

p-mul nthreads dense m n:
    run the dense benchmark on nthreads with powers (m, n)
    mul((1+x+y+z+t)^m, (1+x+y+z+t)^n)
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


static void worker_mul(void * varg)
{
    worker_arg_struct * W = (worker_arg_struct *) varg;
    fmpz_mpoly_mul_threaded(W->G, W->A, W->B, W->ctx, 1);
}

#endif

void profile_mul(
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
    fmpz_mpoly_mul(G, A, B, ctx);
    timeit_stop(timer);
    serial_time = FLINT_MAX(WORD(1), timer->wall);
    flint_printf("    serial time: %wd\n", serial_time);

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
                thread_pool_wake(global_thread_pool, handles[i], 0, worker_mul, worker_args + i);
            }
            else
            {
                worker_mul(worker_args + i);
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

        /* find parallel efficiency */

        fmpz_mpoly_clear(G, ctx);
        fmpz_mpoly_init(G, ctx);
        timeit_start(timer);
        fmpz_mpoly_mul(G, A, B, ctx);
        timeit_stop(timer);
        parallel_time = FLINT_MAX(WORD(1), timer->wall);

        parallel_efficiency = (double)(serial_time)/(double)(parallel_time)/(double)(num_threads);

#if CALCULATE_MACHINE_EFFICIENCY
        flint_printf("parallel %wd time: %wd, efficiency %f (machine %f)\n", num_threads, parallel_time, parallel_efficiency, machine_efficiency);
#else
        flint_printf("parallel %wd time: %wd, efficiency %f\n", num_threads, parallel_time, parallel_efficiency);
#endif
    }

    fmpz_mpoly_clear(G, ctx);
}


int main(int argc, char *argv[])
{
    slong i, m, n, max_threads;
    const slong thread_limit = 64;
    const char * name;

    cpu_affinities = flint_malloc(thread_limit*sizeof(int));
    for (i = 0; i < thread_limit; i++)
        cpu_affinities[i] = i;

    if (argc == 5)
    {
        max_threads = atoi(argv[1]);
        max_threads = FLINT_MIN(max_threads, thread_limit);
        max_threads = FLINT_MAX(max_threads, WORD(1));
        name = argv[2];
        m = atoi(argv[3]);
        n = atoi(argv[4]);
    }
    else
    {
        printf("  usage: p-mul nthreads {dense|sparse} m n\n");
        printf("running: p-mul 4 sparse 12 12\n");
        max_threads = 4;
        name = "sparse";
        m = 12;
        n = 12;
    }

    m = FLINT_MIN(m, WORD(30));
    m = FLINT_MAX(m, WORD(5));
    n = FLINT_MIN(n, WORD(30));
    n = FLINT_MAX(n, WORD(5));

    flint_printf("setting up fmpz_mpoly %s mul ... ", name);

    if (strcmp(name, "dense") == 0)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, A, B;
        const char * vars[] = {"x", "y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(A, ctx);
        fmpz_mpoly_init(B, ctx);

        fmpz_mpoly_set_str_pretty(a, "1 + x + y + z + t", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "1 + x + y + z + t", vars, ctx);
        fmpz_mpoly_pow_ui(A, a, m, ctx);
        fmpz_mpoly_pow_ui(B, b, n, ctx);

        flint_printf("starting dense mul (%wu, %wd):\n", m, n);
        profile_mul(A, B, ctx, max_threads);

        fmpz_mpoly_clear(B, ctx);
        fmpz_mpoly_clear(A, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    else /* "sparse" */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, A, B;
        const char * vars[] = {"x", "y", "z", "t", "u"};

        fmpz_mpoly_ctx_init(ctx, 5, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(A, ctx);
        fmpz_mpoly_init(B, ctx);

        fmpz_mpoly_set_str_pretty(a, "1 + x + y + 2*z^2 + 3*t^3 + 5*u^5", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "1 + u + t + 2*z^2 + 3*y^3 + 5*x^5", vars, ctx);
        fmpz_mpoly_pow_ui(A, a, m, ctx);
        fmpz_mpoly_pow_ui(B, b, n, ctx);

        flint_printf("starting sparse mul (%wu, %wd):\n", m, n);
        profile_mul(A, B, ctx, max_threads);

        fmpz_mpoly_clear(B, ctx);
        fmpz_mpoly_clear(A, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    flint_free(cpu_affinities);
    flint_cleanup_master();
    return 0;
}
