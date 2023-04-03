/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"

int thread_pool_set_size(thread_pool_t T, slong new_size)
{
    thread_pool_entry_struct * D;
    slong i;
    slong old_size;

    new_size = FLINT_MAX(new_size, WORD(0));

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&T->mutex);
#endif
    D = T->tdata;
    old_size = T->length;

    /* check if T is in use */
    for (i = 0; i < old_size; i++)
    {
        if (D[i].available != 1)
        {
#if FLINT_USES_PTHREAD
            pthread_mutex_unlock(&T->mutex);
#endif
            return 0;
        }
    }

    /* destroy all old data */
    for (i = 0; i < old_size; i++)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&D[i].mutex);
#endif
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
    T->tdata = NULL;

    /* create new data */
    if (new_size > 0)
    {
        D = T->tdata
          = (thread_pool_entry_struct *) flint_malloc(new_size
                                           * sizeof(thread_pool_entry_struct));

        for (i = 0; i < new_size; i++)
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

    T->length = new_size;

#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&T->mutex);
#endif
    return 1;
}
