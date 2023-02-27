/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"


void thread_pool_clear(thread_pool_t T)
{
    slong i, size;
    thread_pool_entry_struct * D;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&T->mutex);
#endif
    D = T->tdata;
    size = T->length;
    
    for (i = 0; i < size; i++)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&D[i].mutex);
#endif
	/* all threads should be given back */
        FLINT_ASSERT(D[i].available == 1);
        D[i].exit = 1;
#if FLINT_USES_PTHREAD
	pthread_cond_signal(&D[i].sleep1);
        pthread_mutex_unlock(&D[i].mutex);
        pthread_join(D[i].pth, NULL);
        pthread_cond_destroy(&D[i].sleep2);
        pthread_cond_destroy(&D[i].sleep1);
        pthread_mutex_destroy(&D[i].mutex);
#endif
    }
    if (D != NULL)
    {
        flint_free(D);
    }
#if FLINT_USES_CPUSET && FLINT_USES_PTHREAD
    if (T->original_affinity != NULL) {
        flint_free(T->original_affinity);
        T->original_affinity = NULL;
    }
#endif
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&T->mutex);
    pthread_mutex_destroy(&T->mutex);
#endif
    T->length = -1;
    T->tdata = NULL;
}
