/*
    Copyright 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "fmpz_mpoly.h"
/*
    export LD_LIBRARY_PATH=/tmpbig/schultz/flint2
    likwid-setFrequencies -g performance
    nano nmod_mpoly/profile/p-gcd.c
    make profile MOD=nmod_mpoly
    ./build/nmod_mpoly/profile/p-gcd
*/

slong max_threads;
int * cpu_affinities;

typedef struct _worker_arg_struct
{
    nmod_mpoly_t G;
    const nmod_mpoly_struct * A, * B;
    const nmod_mpoly_ctx_struct * ctx;
} worker_arg_struct;

typedef worker_arg_struct worker_arg_t[1];


static void worker_gcd(void * varg)
{
    worker_arg_struct * W = (worker_arg_struct *) varg;
    nmod_mpoly_gcd_brown_threaded(W->G, W->A, W->B, W->ctx, 1);
}

void profile_gcd(
    const nmod_mpoly_t realG,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const char * name,
    slong m1, slong n1, slong m2, slong n2)
{
    nmod_mpoly_t G;
    timeit_t timer;
    slong num_threads;
    slong serial_time;

    flint_printf("\n******** starting %s ((%wu, %wd), (%wd, %wd)):\n", name, m1, n1, m2, n2);

    flint_set_num_threads(1);
    nmod_mpoly_init(G, ctx);
    timeit_start(timer);
    nmod_mpoly_gcd_brown_threaded(G, A, B, ctx, 1);
    timeit_stop(timer);
    serial_time = FLINT_MAX(WORD(1), timer->wall);
    flint_printf("serial time: %wd\n", serial_time);
    if (!nmod_mpoly_equal(G, realG, ctx))
    {
        printf("gcd wrong\n");
        flint_abort();
    }

    for (num_threads = 2; num_threads <= max_threads; num_threads++)
    {
        thread_pool_handle * handles;
        slong num_workers;
        worker_arg_struct * worker_args;
        slong parallel_time;
        slong i;
        double machine_efficiency, parallel_efficiency;

        flint_set_num_threads(num_threads);
        flint_set_thread_affinity(cpu_affinities, num_threads);

        /* find machine efficiency */

        handles = (thread_pool_handle *) flint_malloc((num_threads - 1)*sizeof(thread_pool_handle));
        num_workers = thread_pool_request(global_thread_pool, handles, num_threads - 1);
        worker_args = (worker_arg_struct *) flint_malloc((num_workers + 1)*sizeof(worker_arg_t));

        timeit_start(timer);
        for (i = 0; i <= num_workers; i++)
        {
            nmod_mpoly_init((worker_args + i)->G, ctx);
            (worker_args + i)->A = A;
            (worker_args + i)->B = B;
            (worker_args + i)->ctx = ctx;
            if (i < num_workers)
            {
                thread_pool_wake(global_thread_pool, handles[i], worker_gcd, worker_args + i);
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
            if (!nmod_mpoly_equal((worker_args + i)->G, realG, ctx))
            {
                printf("gcd wrong\n");
                flint_abort();
            }
            nmod_mpoly_clear((worker_args + i)->G, ctx);

            if (i < num_workers)
            {
                thread_pool_give_back(global_thread_pool, handles[i]);
            }
        }
        flint_free(worker_args);
        flint_free(handles);

        machine_efficiency = (double)(serial_time)/(double)(parallel_time);

        nmod_mpoly_clear(G, ctx);
        nmod_mpoly_init(G, ctx);
        timeit_start(timer);
        nmod_mpoly_gcd_brown_threaded(G, A, B, ctx, MPOLY_DEFAULT_THREAD_LIMIT);
        timeit_stop(timer);
        parallel_time = FLINT_MAX(WORD(1), timer->wall);
        if (!nmod_mpoly_equal(G, realG, ctx))
        {
            printf("gcd wrong\n");
            flint_abort();
        }

        parallel_efficiency = (double)(serial_time)/(double)(parallel_time)/(double)(num_threads);

        flint_printf("parallel %wd time: %wd, efficiency %f (machine %f)\n", num_threads, parallel_time, parallel_efficiency, machine_efficiency);
    }

    nmod_mpoly_clear(G, ctx);
}


int main(int argc, char *argv[])
{
    slong i, m, n, k, l;
    slong step;

    max_threads = argc > 1 ? atoi(argv[1]) : 2;
    max_threads = FLINT_MIN(max_threads, WORD(32));
    max_threads = FLINT_MAX(max_threads, WORD(1));

    flint_printf("setting max_threads = %wd\n", max_threads);

    cpu_affinities = flint_malloc(max_threads*sizeof(int));
    for (i = 0; i < max_threads; i++)
        cpu_affinities[i] = i;

    step = FLINT_MIN(max_threads, 2 + max_threads/7);
    step = FLINT_MAX(step, WORD(3));

    for (m = 0 + 1*step; m <= 0 + 1*step; m += step)
    for (n = 0 + 1*step; n <= 0 + 1*step; n += step)
    for (k = 2 + 0*step; k <= 2 + 1*step; k += step)
    for (l = 2 + 0*step; l <= 2 + 1*step; l += step)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, t, A, B, G;
        const char * vars[] = {"x", "y", "z", "t"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_LEX, n_nextprime(UWORD(1) << (FLINT_BITS - 2), 1));
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_init(A, ctx);
        nmod_mpoly_init(B, ctx);
        nmod_mpoly_init(G, ctx);

        nmod_mpoly_set_str_pretty(a, "1+x+x^2+y^7+z^8+t^40", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "1+x+x^2+y^5+z^9+t^41", vars, ctx);
        nmod_mpoly_pow_ui(A, a, m + k, ctx);
        nmod_mpoly_pow_ui(t, b, n, ctx);
        nmod_mpoly_mul(A, A, t, ctx);
        nmod_mpoly_pow_ui(B, a, m, ctx);
        nmod_mpoly_pow_ui(t, b, n + l, ctx);
        nmod_mpoly_mul(B, B, t, ctx);
        nmod_mpoly_pow_ui(G, a, m, ctx);
        nmod_mpoly_pow_ui(t, b, n, ctx);
        nmod_mpoly_mul(G, G, t, ctx);
        nmod_mpoly_make_monic(G, G, ctx);

        profile_gcd(G, A, B, ctx, "dense gcd", m + k, n, m, n + l);

        nmod_mpoly_clear(G, ctx);
        nmod_mpoly_clear(B, ctx);
        nmod_mpoly_clear(A, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    flint_free(cpu_affinities);

    flint_cleanup();
    return 0;
}
