/*
    Copyright (C) 2013 Fredrik Johansson
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

