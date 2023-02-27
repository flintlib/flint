/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"


void thread_pool_wait(thread_pool_t T, thread_pool_handle i)
{
    thread_pool_entry_struct * D;

    D = T->tdata;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&D[i].mutex);
#endif

    /* should not be trying to wait on a available thread */
    FLINT_ASSERT(D[i].available == 0);

#if FLINT_USES_PTHREAD
    while (D[i].working != 0)
        pthread_cond_wait(&D[i].sleep2, &D[i].mutex);

    pthread_mutex_unlock(&D[i].mutex);
#endif
}
