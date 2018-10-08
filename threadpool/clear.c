/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "threadpool.h"


void threadpool_clear(threadpool_t T)
{
    slong i;
    tpentry_struct * D = T->tdata;
    slong l = T->length;

    for (i = 0; i < l; i++)
    {
        FLINT_ASSERT(D[i].available == 1); /* all threads should be given back */
        pthread_mutex_lock(&D[i].mutex);
        D[i].exit = 1;
        pthread_cond_signal(&D[i].sleep1);
        pthread_mutex_unlock(&D[i].mutex);
        pthread_join(D[i].pth, NULL);
        pthread_cond_destroy(&D[i].sleep2);
        pthread_cond_destroy(&D[i].sleep1);
        pthread_mutex_destroy(&D[i].mutex);
    }
    flint_free(D);
    pthread_mutex_destroy(&T->mutex);
    T->length = -1;
    T->tdata = NULL;
}
