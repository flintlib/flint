/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <pthread.h>
#include "flint.h"

typedef struct
{
    pthread_t pth;
    pthread_mutex_t mutex;
    pthread_cond_t sleep1;
    pthread_cond_t sleep2;
    volatile int idx;
    volatile int available;
    void (*fxn)(void *); /* <- TODO:                                        */
    void * fxnarg;       /* <-   these should be marked volatile somehow ?? */
    volatile int working;
    volatile int exit;
} tpentry_struct;

typedef tpentry_struct tpentry_t[1];


typedef struct
{
    pthread_mutex_t mutex;
    tpentry_struct * tdata;
    slong length;
} threadpool_struct;

typedef threadpool_struct threadpool_t[1];

typedef int threadpool_threadhandle;
