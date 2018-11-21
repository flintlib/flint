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
    void (* fxn)(void *);
    void * fxnarg;
    volatile int working;
    volatile int exit;
} thread_pool_entry_struct;

typedef thread_pool_entry_struct thread_pool_entry_t[1];

typedef struct
{
    pthread_mutex_t mutex;
    thread_pool_entry_struct * tdata;
    slong length;
} thread_pool_struct;

typedef thread_pool_struct thread_pool_t[1];

typedef int thread_pool_handle;

extern thread_pool_t global_thread_pool;
extern int global_thread_pool_initialized;

FLINT_DLL void * thread_pool_idle_loop(void * varg);

FLINT_DLL void thread_pool_init(thread_pool_t T, slong l);

FLINT_DLL slong thread_pool_get_size(thread_pool_t T);

FLINT_DLL int thread_pool_set_size(thread_pool_t T, slong new_size);

FLINT_DLL slong thread_pool_request(thread_pool_t T,
                                    thread_pool_handle * out, slong requested);

FLINT_DLL void thread_pool_wake(thread_pool_t T, thread_pool_handle i,
                                                   void (*f)(void*), void * a);

FLINT_DLL void thread_pool_wait(thread_pool_t T, thread_pool_handle i);

FLINT_DLL void thread_pool_give_back(thread_pool_t T, thread_pool_handle i);

FLINT_DLL void thread_pool_clear(thread_pool_t T);
