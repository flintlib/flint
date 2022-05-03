/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include "flint.h"

#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
#if FLINT_USES_PTHREAD
    pthread_t pth;
    pthread_mutex_t mutex;
    pthread_cond_t sleep1;
    pthread_cond_t sleep2;
#endif
    volatile int idx;
    volatile int available;
    volatile int max_workers;
    void (* fxn)(void *);
    void * fxnarg;
    volatile int working;
    volatile int exit;
} thread_pool_entry_struct;

typedef thread_pool_entry_struct thread_pool_entry_t[1];

typedef struct
{
#if FLINT_USES_CPUSET && FLINT_USES_PTHREAD
    void * original_affinity;
#endif
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    thread_pool_entry_struct * tdata;
    slong length;
} thread_pool_struct;

typedef thread_pool_struct thread_pool_t[1];

typedef int thread_pool_handle;

FLINT_DLL extern thread_pool_t global_thread_pool;
FLINT_DLL extern int global_thread_pool_initialized;

FLINT_DLL void * thread_pool_idle_loop(void * varg);

FLINT_DLL void thread_pool_init(thread_pool_t T, slong l);

FLINT_DLL int thread_pool_set_affinity(thread_pool_t T,
                                                     int * cpus, slong length);

FLINT_DLL int thread_pool_restore_affinity(thread_pool_t T);

FLINT_DLL slong thread_pool_get_size(thread_pool_t T);

FLINT_DLL int thread_pool_set_size(thread_pool_t T, slong new_size);

FLINT_DLL slong thread_pool_request(thread_pool_t T,
                                    thread_pool_handle * out, slong requested);

FLINT_DLL void thread_pool_wake(thread_pool_t T, thread_pool_handle i,
                                  int max_workers, void (*f)(void*), void * a);

FLINT_DLL void thread_pool_wait(thread_pool_t T, thread_pool_handle i);

FLINT_DLL void thread_pool_give_back(thread_pool_t T, thread_pool_handle i);

FLINT_DLL void thread_pool_clear(thread_pool_t T);

/* misc internal helpers *****************************************************/

FLINT_DLL void _thread_pool_distribute_work_2(slong start, slong stop,
                                    slong * Astart, slong * Astop, slong Alen,
                                    slong * Bstart, slong * Bstop, slong Blen);

FLINT_DLL ulong _thread_pool_find_work_2(ulong a, ulong alpha,
                                      ulong b, ulong beta, ulong yn, ulong yd);

#ifdef __cplusplus
}
#endif

#endif
