/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "threadpool.h"


void threadpool_wake(threadpool_t T, threadpool_threadhandle i,
                                                    void (*f)(void*), void * a)
{
    tpentry_struct * D;

    pthread_mutex_lock(&T->mutex);

    FLINT_ASSERT(i < T->length);

    D = T->tdata;

    FLINT_ASSERT(D[i].available == 0); /* should not be trying to wake an available thread */

    pthread_mutex_lock(&D[i].mutex);
    D[i].working = 1;
    D[i].fxn = f;
    D[i].fxnarg = a;
    pthread_cond_signal(&D[i].sleep1);
    pthread_mutex_unlock(&D[i].mutex);

    pthread_mutex_unlock(&T->mutex);
}
