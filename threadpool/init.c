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

slong threadpool_size(threadpool_t T)
{
    slong ret;
    pthread_mutex_lock(&T->mutex);
    ret = T->length;
    pthread_mutex_unlock(&T->mutex);
    return ret;
}

slong threadpool_request(threadpool_t T, threadpool_threadhandle * out, slong requested)
{
    slong i, ret = 0;
    tpentry_struct * D;

    if (requested <= 0)
        return 0;

    pthread_mutex_lock(&T->mutex);

    D = T->tdata;
    if (T->length > 0)
    {
        for (i = 0; i < T->length; i++)
        {
            if (D[i].available == 1)
            {
                D[i].available = 0;
                out[ret] = i;
                ret++;
                if (ret >= requested)
                    break;
            }
        }
    }

    pthread_mutex_unlock(&T->mutex);

    return ret;
}

void threadpool_wake(threadpool_t T, threadpool_threadhandle i, void (*f)(void*), void * a)
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

void threadpool_wait(threadpool_t T, threadpool_threadhandle i)
{
    tpentry_struct * D;

    D = T->tdata;

    FLINT_ASSERT(D[i].available == 0); /* should not be trying to wait on a available thread */

    pthread_mutex_lock(&D[i].mutex);
    while (D[i].working != 0)
        pthread_cond_wait(&D[i].sleep2, &D[i].mutex);        
    pthread_mutex_unlock(&D[i].mutex);

    pthread_mutex_lock(&T->mutex);
    D[i].available = 0;
    pthread_mutex_unlock(&T->mutex);
}

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

void threadpool_clear(threadpool_t T)
{
    slong i;
    tpentry_struct * D = T->tdata;
    slong l = T->length;

    for (i = 0; i < l; i++)
    {
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
