/*
    Copyright (C) 2013, 2022 Fredrik Johansson
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "thread_pool.h"
#include "thread_support.h"

/* Automatically initialised to zero when threads are started */
FLINT_TLS_PREFIX int _flint_num_workers = 0;

int flint_get_num_threads()
{
    return _flint_num_workers + 1;
}

void flint_set_num_threads(int num_threads)
{
#if !FLINT_USES_PTHREAD
    num_threads = 1;
#endif
    _flint_num_workers = num_threads - 1;
    if (global_thread_pool_initialized)
    {
        if (!thread_pool_set_size(global_thread_pool, num_threads - 1))
        {
            flint_throw(FLINT_ERROR,
               "flint_set_num_threads called while global thread pool in use");
        }
    }
    else
    {
        thread_pool_init(global_thread_pool, num_threads - 1);
        global_thread_pool_initialized = 1;
    }
}

void _flint_set_num_workers(int num_workers)
{
    _flint_num_workers = num_workers;
}

int flint_set_num_workers(int num_workers)
{
    int old_num_workers;

#if !FLINT_USES_PTHREAD
    num_workers = 0;
#endif

    old_num_workers = _flint_num_workers;
    
    _flint_num_workers = FLINT_MIN(_flint_num_workers, num_workers);

    return old_num_workers;
}

void flint_reset_num_workers(int num_workers)
{
    _flint_num_workers = num_workers;
}

/* return zero for success, nonzero for error */
int flint_set_thread_affinity(int * cpus, slong length)
{
    if (!global_thread_pool_initialized)
        return 1;

    return thread_pool_set_affinity(global_thread_pool, cpus, length);
}

/* return zero for success, nonzero for error */
int flint_restore_thread_affinity()
{
    if (!global_thread_pool_initialized)
        return 1;

    return thread_pool_restore_affinity(global_thread_pool);
}

slong flint_request_threads(thread_pool_handle ** handles, slong thread_limit)
{
    slong num_handles = 0;
    slong num_threads = flint_get_num_threads();
    
    thread_limit = FLINT_MIN(thread_limit, num_threads);

    *handles = NULL;

    if (global_thread_pool_initialized && thread_limit > 1)
    {
        slong max_num_handles;
        max_num_handles = thread_pool_get_size(global_thread_pool);
        max_num_handles = FLINT_MIN(thread_limit - 1, max_num_handles);
        if (max_num_handles > 0)
        {
            *handles = (thread_pool_handle *) flint_malloc(
                                   max_num_handles*sizeof(thread_pool_handle));
            num_handles = thread_pool_request(global_thread_pool,
                                                     *handles, max_num_handles);
        }
    }

    return num_handles;
}

void flint_give_back_threads(thread_pool_handle * handles, slong num_handles)
{
    slong i;

    for (i = 0; i < num_handles; i++)
        thread_pool_give_back(global_thread_pool, handles[i]);

    if (handles)
        flint_free(handles);
}

typedef struct
{
    do_func_t f;
    void * args;
    slong a;
    slong b;
    slong step;
}
work_chunk_t;

static void
worker(void * _work)
{
    work_chunk_t work = *((work_chunk_t *) _work);
    slong i;

    for (i = work.a; i < work.b; i += work.step)
        work.f(i, work.args);
}

void flint_parallel_do(do_func_t f, void * args, slong n, int thread_limit, int flags)
{
    slong i;

    if (thread_limit <= 0)
        thread_limit = flint_get_num_threads();

    thread_limit = FLINT_MIN(thread_limit, n);

    if (thread_limit <= 1)
    {
        for (i = 0; i < n; i++)
            f(i, args);
    }
    else
    {
        slong i, num_threads, num_workers;
        thread_pool_handle * handles;

        num_workers = flint_request_threads(&handles, thread_limit);
        num_threads = num_workers + 1;

        if (flags & FLINT_PARALLEL_VERBOSE)
            flint_printf("parallel_do with num_threads = %wd\n", num_threads);

        if (num_workers < 1)
        {
            flint_give_back_threads(handles, num_workers);

            for (i = 0; i < n; i++)
                f(i, args);
        }
        else
        {
            work_chunk_t * work;
            slong chunk_size;
            TMP_INIT;
            TMP_START;

            work = TMP_ALLOC(num_threads * sizeof(work_chunk_t));

            if (flags & FLINT_PARALLEL_STRIDED)
            {
                for (i = 0; i < num_threads; i++)
                {
                    work[i].f = f;
                    work[i].args = args;
                    work[i].a = i;
                    work[i].b = n;
                    work[i].step = num_threads;
                }
            }
            else
            {
                chunk_size = (n + num_threads - 1) / num_threads;

                for (i = 0; i < num_threads; i++)
                {
                    work[i].f = f;
                    work[i].args = args;
                    work[i].a = i * chunk_size;
                    work[i].b = FLINT_MIN((i + 1) * chunk_size, n);
                    work[i].step = 1;
                }
            }

            if (flags & FLINT_PARALLEL_VERBOSE)
            {
                for (i = 0; i < num_threads; i++)
                {
                    flint_printf("thread #%wd allocated a = %wd, b = %wd, step = %wd\n", i, work[i].a, work[i].b, work[i].step);
                }
            }

            for (i = 0; i < num_workers; i++)
                thread_pool_wake(global_thread_pool, handles[i], 0, worker, &work[i]);

            worker(&work[num_workers]);

            for (i = 0; i < num_workers; i++)
                thread_pool_wait(global_thread_pool, handles[i]);

            flint_give_back_threads(handles, num_workers);
            TMP_END;
        }
    }
}

typedef struct
{
    void * res;
    bsplit_basecase_func_t basecase;
    bsplit_merge_func_t merge;
    size_t sizeof_res;
    bsplit_init_func_t init;
    bsplit_clear_func_t clear;
    void * args;
    slong a;
    slong b;
    slong basecase_cutoff;
    slong thread_limit;
    int flags;
}
flint_parallel_binary_splitting_t;

static void
_bsplit_worker(void * _args)
{
    flint_parallel_binary_splitting_t * args = (flint_parallel_binary_splitting_t *) _args;

    flint_parallel_binary_splitting(args->res, args->basecase, args->merge, args->sizeof_res, args->init, args->clear, args->args, args->a, args->b, args->basecase_cutoff, args->thread_limit, args->flags);
}

void
flint_parallel_binary_splitting(void * res, bsplit_basecase_func_t basecase, bsplit_merge_func_t merge,
    size_t sizeof_res, bsplit_init_func_t init, bsplit_clear_func_t clear, void * args, slong a, slong b, slong basecase_cutoff, int thread_limit, int flags)
{
    basecase_cutoff = FLINT_MAX(basecase_cutoff, 1);

    if (b - a <= basecase_cutoff)
    {
        basecase(res, a, b, args);
    }
    else
    {
        void * left, * right;
        slong m = a + (b - a) / 2;
        slong nw;
        slong nw_save;
        slong nt;
        thread_pool_handle * threads;
        TMP_INIT;

        TMP_START;

        if (flags & FLINT_PARALLEL_BSPLIT_LEFT_INPLACE)
        {
            left = res;
            right = TMP_ALLOC(sizeof_res);

            init(right, args);
        }
        else
        {
            left = TMP_ALLOC(2 * sizeof_res);
            right = (void *) (((char *) left) + sizeof_res);

            init(left, args);
            init(right, args);
        }

        if (thread_limit <= 0)
            thread_limit = flint_get_num_threads();

        nt = thread_limit;
        nw = flint_request_threads(&threads, FLINT_MIN(nt, 2)); /* request one extra worker */

        if (nw == 0)
        {
            flint_parallel_binary_splitting(left, basecase, merge, sizeof_res, init, clear, args, a, m, basecase_cutoff, thread_limit, flags);
            flint_parallel_binary_splitting(right, basecase, merge, sizeof_res, init, clear, args, m, b, basecase_cutoff, thread_limit, flags);
        }
        else
        {
            flint_parallel_binary_splitting_t right_args;

            FLINT_ASSERT(nt >= 2);

            nw_save = flint_set_num_workers(nt - nt / 2 - 1);

            right_args.res = right;
            right_args.basecase = basecase;
            right_args.merge = merge;
            right_args.sizeof_res = sizeof_res;
            right_args.init = init;
            right_args.clear = clear;
            right_args.args = args;
            right_args.a = m;
            right_args.b = b;
            right_args.basecase_cutoff = basecase_cutoff;
            right_args.thread_limit = thread_limit;
            right_args.flags = flags;

            thread_pool_wake(global_thread_pool, threads[0], nt / 2 - 1, _bsplit_worker, &right_args);

            flint_parallel_binary_splitting(left, basecase, merge, sizeof_res, init, clear, args, a, m, basecase_cutoff, thread_limit, flags);

            flint_reset_num_workers(nw_save);
            thread_pool_wait(global_thread_pool, threads[0]);
        }

        flint_give_back_threads(threads, nw);

        merge(res, left, right, args);

        if (flags & FLINT_PARALLEL_BSPLIT_LEFT_INPLACE)
        {
            clear(right, args);
        }
        else
        {
            clear(left, args);
            clear(right, args);
        }

        TMP_END;
    }
}
