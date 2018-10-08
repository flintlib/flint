/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "threadpool.h"


void * threadpool_idle_loop(void * varg)
{
    tpentry_struct * arg = (tpentry_struct *) varg;

    goto threadpool_Lock;

threadpool_DoWork:

    arg->fxn(arg->fxnarg);

threadpool_Lock:

    pthread_mutex_lock(&arg->mutex);
    arg->working = 0;

threadpool_CheckExit:

    if (arg->exit != 0)
        goto threadpool_Unlock;

    pthread_cond_signal(&arg->sleep2);
    pthread_cond_wait(&arg->sleep1, &arg->mutex);

    if (arg->working == 0)
        goto threadpool_CheckExit;

threadpool_Unlock:

    pthread_mutex_unlock(&arg->mutex);

    if (arg->exit == 0)
        goto threadpool_DoWork;

    flint_cleanup();

    return NULL;
}


void threadpool_init(threadpool_t T, slong l)
{
    slong i;
    tpentry_struct * D;
    FLINT_ASSERT(l >= 0);

    pthread_mutex_init(&T->mutex, NULL);
    T->length = l;

    if (l == 0)
    {
        T->tdata = NULL;
        return;
    }

    D = (tpentry_struct *) flint_malloc(l * sizeof(tpentry_struct));
    T->tdata = D;

    for (i = 0; i < l; i++)
    {
        pthread_mutex_init(&D[i].mutex, NULL);
        pthread_cond_init(&D[i].sleep1, NULL);
        pthread_cond_init(&D[i].sleep2, NULL);
        D[i].idx = i;
        D[i].available = 1;
        D[i].fxn = NULL;
        D[i].fxnarg = NULL;
        D[i].working = -1;
        D[i].exit = 0;
        pthread_mutex_lock(&D[i].mutex);
        pthread_create(&D[i].pth, NULL, threadpool_idle_loop, &D[i]);
        while (D[i].working != 0)
            pthread_cond_wait(&D[i].sleep2, &D[i].mutex);
        pthread_mutex_unlock(&D[i].mutex);
    }
}
