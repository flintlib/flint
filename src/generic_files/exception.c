/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "flint.h"

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
#include <pthread.h>

static pthread_once_t abort_func_init = PTHREAD_ONCE_INIT;
pthread_mutex_t abort_func_lock;

void __flint_set_abort_init(void)
{
   pthread_mutex_init(&abort_func_lock, NULL);
}
#endif

#ifdef __GNUC__
__attribute__((noreturn)) void (*abort_func)(void) = abort;
#else
void (*abort_func)(void) = abort;
#endif

void flint_set_abort(void (*func)(void))
{
#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&abort_func_init, __flint_set_abort_init);
    pthread_mutex_lock(&abort_func_lock);
#endif

#ifdef __GNUC__
# pragma GCC diagnostic push
# if defined(__clang__)
#  pragma GCC diagnostic ignored "-Wincompatible-function-pointer-types"
# elif defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wdiscarded-qualifiers"
# endif
#endif
    abort_func = func;
#ifdef __GNUC__
# pragma GCC diagnostic pop
#endif

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&abort_func_lock);
#endif
}

FLINT_NORETURN void flint_abort(void)
{
    fflush(stdout);
    fflush(stderr);
    (*abort_func)();
}

void flint_throw(flint_err_t exc, const char * msg, ...)
{
    va_list ap;

    va_start(ap, msg);

    if (exc != FLINT_TEST_FAIL)
    {
        printf("Flint exception (");
        switch (exc)
        {
            case FLINT_ERROR:
                printf("General error");
                break;
            case FLINT_OVERFLOW:
                printf("Overflow");
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
        printf("):\n    ");
    }
    else
    {
        printf("FAIL!\n\n");
    }

    flint_vprintf(msg, ap);
    va_end(ap);

    flint_abort();
}
