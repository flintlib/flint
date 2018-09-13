/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "threadpool.h"


void threadpool_giveback(threadpool_t T, threadpool_threadhandle i)
{
    tpentry_struct * D;

    D = T->tdata;

    FLINT_ASSERT(D[i].available == 0); /* should not be trying to giveback an available thread */

    pthread_mutex_lock(&D[i].mutex);
    FLINT_ASSERT(D[i].working == 0);
    pthread_mutex_unlock(&D[i].mutex);

    pthread_mutex_lock(&T->mutex);
    D[i].available = 1;    
    pthread_mutex_unlock(&T->mutex);
}
