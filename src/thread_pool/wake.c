/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"


void thread_pool_wake(thread_pool_t T, thread_pool_handle i,
                                   int max_workers, void (*f)(void*), void * a)
{
    thread_pool_entry_struct * D;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&T->mutex);
#endif

    FLINT_ASSERT(i < T->length);

    D = T->tdata;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&D[i].mutex);
#endif

    /* should not be trying to wake an available thread */
    FLINT_ASSERT(D[i].available == 0);

    D[i].max_workers = max_workers;
    D[i].working = 1;
    D[i].fxn = f;
    D[i].fxnarg = a;
#if FLINT_USES_PTHREAD
    pthread_cond_signal(&D[i].sleep1);

    pthread_mutex_unlock(&D[i].mutex);

    pthread_mutex_unlock(&T->mutex);
#endif
}
