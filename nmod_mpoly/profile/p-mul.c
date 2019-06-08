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
    nano nmod_mpoly/profile/p-mul.c
    make profile MOD=nmod_mpoly
    ./build/nmod_mpoly/profile/p-mul
*/

slong max_threads;
int * cpu_affinities;

typedef struct _worker_arg_struct
{
    nmod_mpoly_t P;
    const nmod_mpoly_struct * A, * B;
    const nmod_mpoly_ctx_struct * ctx;
} worker_arg_struct;

typedef worker_arg_struct worker_arg_t[1];


static void worker_mul(void * varg)
{
    worker_arg_struct * W = (worker_arg_struct *) varg;
    nmod_mpoly_mul_threaded(W->P, W->A, W->B, W->ctx, 1);
}

void profile_mul(
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const char * name,
    slong m1, slong n1)
{
    nmod_mpoly_t P;
    timeit_t timer;
    slong num_threads;
    slong serial_time;

    flint_printf("\n******** starting %s (%wu, %wd):\n", name, m1, n1);

    flint_set_num_threads(1);
    nmod_mpoly_init(P, ctx);
    timeit_start(timer);
    nmod_mpoly_mul(P, A, B, ctx);
    timeit_stop(timer);
    serial_time = FLINT_MAX(WORD(1), timer->wall);
    flint_printf("    serial time: %wd\n", serial_time);

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

        handles = (thread_pool_handle *) flint_malloc((num_threads - 1)*sizeof(thread_pool_handle));
        num_workers = thread_pool_request(global_thread_pool, handles, num_threads - 1);
        worker_args = (worker_arg_struct *) flint_malloc((num_workers + 1)*sizeof(worker_arg_t));

        timeit_start(timer);
        for (i = 0; i <= num_workers; i++)
        {
            nmod_mpoly_init((worker_args + i)->P, ctx);
            (worker_args + i)->A = A;
            (worker_args + i)->B = B;
            (worker_args + i)->ctx = ctx;
            if (i < num_workers)
            {
                thread_pool_wake(global_thread_pool, handles[i], worker_mul, worker_args + i);
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
            nmod_mpoly_clear((worker_args + i)->P, ctx);

            if (i < num_workers)
            {
                thread_pool_give_back(global_thread_pool, handles[i]);
            }
        }
        flint_free(worker_args);
        flint_free(handles);

        machine_efficiency = (double)(serial_time)/(double)(parallel_time);

        /* find parallel efficiency */

        nmod_mpoly_clear(P, ctx);
        nmod_mpoly_init(P, ctx);
        timeit_start(timer);
        nmod_mpoly_mul(P, A, B, ctx);
        timeit_stop(timer);
        parallel_time = FLINT_MAX(WORD(1), timer->wall);

        parallel_efficiency = (double)(serial_time)/(double)(parallel_time)/(double)(num_threads);

        flint_printf("parallel %wd time: %wd, efficiency %f (machine %f)\n", num_threads, parallel_time, parallel_efficiency, machine_efficiency);
    }

    nmod_mpoly_clear(P, ctx);
}


int main(int argc, char *argv[])
{
    slong i, m, n;

    max_threads = argc > 1 ? atoi(argv[1]) : 2;
    max_threads = FLINT_MIN(max_threads, WORD(32));
    max_threads = FLINT_MAX(max_threads, WORD(1));

    flint_printf("setting max_threads = %wd\n", max_threads);

    cpu_affinities = flint_malloc(max_threads*sizeof(int));
    for (i = 0; i < max_threads; i++)
        cpu_affinities[i] = i;

    for (m = 6 + max_threads/4; m <= 12 + max_threads/4; m += 3)
    for (n = 6 + max_threads/4; n <= 12 + max_threads/4; n += 3)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, A, B;
        const char * vars[] = {"x", "y", "z", "t", "u"};

        nmod_mpoly_ctx_init(ctx, 5, ORD_LEX, n_nextprime(UWORD(1) << (FLINT_BITS - 2), 1));
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(A, ctx);
        nmod_mpoly_init(B, ctx);

        nmod_mpoly_set_str_pretty(a, "1+x+1*y^2+1*z^3+1*t^4+1*u^5", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "1+u+1*t^2+1*z^3+1*y^4+1*x^5", vars, ctx);
        nmod_mpoly_pow_ui(A, a, m, ctx);
        nmod_mpoly_pow_ui(B, b, n, ctx);

        profile_mul(A, B, ctx, "sparse mul", m, n);

        nmod_mpoly_clear(B, ctx);
        nmod_mpoly_clear(A, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    flint_free(cpu_affinities);

    flint_cleanup();
    return 0;
}
