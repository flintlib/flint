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

    flint_printf("Flint exception (");

    switch (exc)
    {
        case FLINT_ERROR:
            flint_printf("General error");
            break;
        case FLINT_IMPINV:
            flint_printf("Impossible inverse");
            break;
        case FLINT_DOMERR:
            flint_printf("Domain error");
            break;
        case FLINT_DIVZERO:
            flint_printf("Divide by zero");
            break;
        case FLINT_INEXACT:
            flint_printf("Inexact");
            break;
        default:
            flint_printf("Unknown exception");
     }

     printf("):\n    ");

     flint_vprintf(msg, ap);
     va_end(ap);

     flint_abort();
}

/**** overridable exception function *****************************************/

FLINT_NORETURN static void
def_exception_func(int exc, const char * msg, char * extra)
{
    flint_printf("flint exception (");

    if (exc == FLINT_ERROR)
        flint_printf("General error");
    else if (exc == FLINT_IMPINV)
        flint_printf("Impossible inverse");
    else if (exc == FLINT_DOMERR)
        flint_printf("Domain error");
    else if (exc == FLINT_DIVZERO)
        flint_printf("Divide by zero");
    else if (exc == FLINT_INEXACT)
        flint_printf("Inexact");
    else
        flint_printf("Unknown");

    flint_printf("): ");
    flint_printf(msg);

    if (extra != NULL)
    {
        flint_printf(extra);
        flint_free(extra);
    }

    flint_printf("\n");
    fflush(stdout);

    flint_abort();
}

FLINT_NORETURN void (*exception_func)(int, const char *, char *) = def_exception_func;

FLINT_NORETURN void
flint_exception(flint_err_t exc, const char * msg, char * extra)
{
    (*exception_func)((int)(exc), msg, extra);
}

void flint_set_exception(FLINT_NORETURN void (*func)(int, const char *, char *))
{
    exception_func = func;
}

