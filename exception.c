/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx
#include <stdio.h>
#include <stdlib.h>
#undef ulong
#include "flint.h"

#if FLINT_REENTRANT && !FLINT_USES_TLS
#include <pthread.h>

static pthread_once_t abort_func_init = PTHREAD_ONCE_INIT;
pthread_mutex_t abort_func_lock;

void __flint_set_abort_init()
{
    pthread_mutex_init(&abort_func_lock, NULL);
}
#endif

FLINT_NORETURN void (*abort_func)(void) = abort;

void flint_set_abort(FLINT_NORETURN void (*func)(void))
{
#if FLINT_REENTRANT && !FLINT_USES_TLS
    pthread_once(&abort_func_init, __flint_set_abort_init);
    pthread_mutex_lock(&abort_func_lock);
#endif

    abort_func = func;

#if FLINT_REENTRANT && !FLINT_USES_TLS
    pthread_mutex_unlock(&abort_func_lock);
#endif
}

FLINT_NORETURN void flint_abort()
{
    (*abort_func)();
}

void flint_throw(flint_err_t exc, const char * msg, ...)
{
    printf("Flint exception (");
    switch (exc)
    {
        case FLINT_ERROR:
            printf("General error");
            break;
        case FLINT_ALLOC:
            printf("Allocation error");
            break;
        case FLINT_MEMMGR:
            printf("Memory manager error");
            break;
        case FLINT_IMPINV:
            printf("Impossible inverse");
            break;
        case FLINT_DOMERR:
            printf("Domain error");
            break;
        case FLINT_DIVZERO:
            printf("Divide by zero");
            break;
        case FLINT_EXPOF:
            printf("Exponent overflow");
            break;
        case FLINT_INEXACT:
            printf("Inexact");
            break;
        default:
            printf("Unknown exception");
    }
    printf("):  ");

    printf(msg);
    fflush(stdout);

    flint_abort();
}
