/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
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
    va_list ap;

    va_start(ap, msg);

    switch (exc)
    {
        case FLINT_ERROR:
            printf("Flint exception (General error):\n    ");
            break;
        case FLINT_IMPINV:
            printf("Flint exception (Impossible inverse):\n    ");
            break;
        case FLINT_DOMERR:
            printf("Flint exception (Domain error):\n    ");
            break;
        case FLINT_DIVZERO:
            printf("Flint exception (Divide by zero):\n    ");
            break;
        case FLINT_INEXACT:
            printf("Flint exception (Inexact):\n    ");
            break;
        default:
            printf("Flint exception (Unknown exception):\n    ");
     }

     flint_vprintf(msg, ap);
     va_end(ap);

     flint_abort();
}
