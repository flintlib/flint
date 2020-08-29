/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"


void thread_pool_give_back(thread_pool_t T, thread_pool_handle i)
{
    thread_pool_entry_struct * D;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&T->mutex);
#endif
    D = T->tdata;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&D[i].mutex);
#endif

    /* thread we are giving back should not be available nor working */
    FLINT_ASSERT(D[i].available == 0);
    FLINT_ASSERT(D[i].working == 0);

    D[i].available = 1;

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&D[i].mutex);

    pthread_mutex_unlock(&T->mutex);
#endif
}
