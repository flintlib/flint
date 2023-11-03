/*
    Copyright (C) 2018-2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#include <sched.h>

#include "thread_pool.h"

#if FLINT_USES_PTHREAD && defined(HAVE_PTHREAD_NP_H)
# include <pthread_np.h>
#endif

thread_pool_t global_thread_pool;
int global_thread_pool_initialized = 0;


void * thread_pool_idle_loop(void * varg)
{
    thread_pool_entry_struct * arg = (thread_pool_entry_struct *) varg;

    goto thread_pool_Lock;

thread_pool_DoWork:

    _flint_set_num_workers(arg->max_workers);
    arg->fxn(arg->fxnarg);

thread_pool_Lock:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&arg->mutex);
#endif
    arg->working = 0;

thread_pool_CheckExit:

    if (arg->exit != 0)
        goto thread_pool_Unlock;

#if FLINT_USES_PTHREAD
    pthread_cond_signal(&arg->sleep2);
    pthread_cond_wait(&arg->sleep1, &arg->mutex);
#endif

    if (arg->working == 0)
        goto thread_pool_CheckExit;

thread_pool_Unlock:

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&arg->mutex);
#endif

    if (arg->exit == 0)
        goto thread_pool_DoWork;

    flint_cleanup();

    return NULL;
}


void thread_pool_init(thread_pool_t T, slong size)
{
    slong i;
    thread_pool_entry_struct * D;
    size = FLINT_MAX(size, WORD(0));

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&T->mutex, NULL);
#endif
    T->length = size;

#if FLINT_USES_CPUSET && FLINT_USES_PTHREAD
    T->original_affinity = flint_malloc(sizeof(cpu_set_t));
    if (0 != pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t),
                                        (cpu_set_t *)T->original_affinity))
    {
        CPU_ZERO((cpu_set_t *)T->original_affinity);
    }
#endif

    if (size == 0)
    {
        T->tdata = NULL;
        return;
    }

    D = (thread_pool_entry_struct *) flint_malloc(
                                      size * sizeof(thread_pool_entry_struct));
    T->tdata = D;

    for (i = 0; i < size; i++)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_init(&D[i].mutex, NULL);
        pthread_cond_init(&D[i].sleep1, NULL);
        pthread_cond_init(&D[i].sleep2, NULL);
#endif
	D[i].idx = i;
        D[i].available = 1;
        D[i].fxn = NULL;
        D[i].fxnarg = NULL;
        D[i].working = -1;
	D[i].max_workers = 0;
        D[i].exit = 0;
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&D[i].mutex);
        pthread_create(&D[i].pth, NULL, thread_pool_idle_loop, &D[i]);
	while (D[i].working != 0)
            pthread_cond_wait(&D[i].sleep2, &D[i].mutex);
        pthread_mutex_unlock(&D[i].mutex);
#endif
    }
}
