/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "flint.h"

#if HAVE_GC
#include "gc.h"
#endif

#if FLINT_REENTRANT && !HAVE_TLS
#include <pthread.h>

static pthread_once_t register_initialised = PTHREAD_ONCE_INIT;
pthread_mutex_t register_lock;
#endif

static __inline__ void * __default_flint_malloc(size_t size)
{
#if HAVE_GC
    return GC_malloc(size);
#else
    return malloc(size);
#endif
}

static __inline__ void * __default_flint_calloc(size_t num, size_t size)
{
#if HAVE_GC
    return GC_malloc(num*size);
#else
    return calloc(num, size);
#endif
}

static __inline__ void * __default_flint_realloc(void * ptr, size_t size)
{
#if HAVE_GC
    return GC_realloc(ptr, size);
#else
    return realloc(ptr, size);
#endif
}

static __inline__ void __default_flint_free(void * ptr)
{
#if !HAVE_GC
    free(ptr);
#endif
}

static void * (*__flint_malloc) (size_t size) = __default_flint_malloc;
static void * (*__flint_realloc) (void * ptr, size_t size) = __default_flint_realloc;
static void * (*__flint_calloc) (size_t num, size_t size) = __default_flint_calloc;
static void (*__flint_free) (void * ptr) = __default_flint_free;

static void flint_memory_error()
{
    flint_printf("Exception (FLINT memory_manager). Unable to allocate memory.\n");
    abort();
}

void flint_set_memory_functions (
    void * (*flint_malloc) (size_t size),
    void * (*flint_realloc) (void * ptr, size_t size),
    void * (*flint_calloc) (size_t num, size_t size),
    void (*flint_free) (void * ptr))
{
  if (flint_malloc == 0)
    flint_malloc = __default_flint_malloc;
  if (flint_realloc == 0)
    flint_realloc = __default_flint_realloc;
  if (flint_calloc == 0)
    flint_calloc = __default_flint_calloc;
  if (flint_free == 0)
    flint_free = __default_flint_free;

  __flint_malloc = flint_malloc;
  __flint_realloc = flint_realloc;
  __flint_calloc = flint_calloc;
  __flint_free = flint_free;
}

void * flint_malloc(size_t size)
{
    void * ptr;

    ptr = __flint_malloc(size);

    if (ptr == NULL)
        flint_memory_error();

    return ptr;
}

void * flint_realloc(void * ptr, size_t size)
{
    void * ptr2;

    ptr2 = __flint_realloc(ptr, size);

    if (ptr2 == NULL)
        flint_memory_error();

    return ptr2;
}

void * flint_calloc(size_t num, size_t size)
{
    void * ptr;

    ptr = __flint_calloc(num, size);

    if (ptr == NULL)
        flint_memory_error();

    return ptr;
}

void flint_free(void * ptr)
{
    __flint_free(ptr);
}


FLINT_TLS_PREFIX size_t flint_num_cleanup_functions = 0;

FLINT_TLS_PREFIX flint_cleanup_function_t * flint_cleanup_functions = NULL;

#if FLINT_REENTRANT && !HAVE_TLS
void register_init()
{
   pthread_mutex_init(&register_lock, NULL);
}
#endif

void flint_register_cleanup_function(flint_cleanup_function_t cleanup_function)
{
#if FLINT_REENTRANT && !HAVE_TLS
    pthread_once(&register_initialised, register_init);
    pthread_mutex_lock(&register_lock);
#endif

    flint_cleanup_functions = flint_realloc(flint_cleanup_functions,
        (flint_num_cleanup_functions + 1) * sizeof(flint_cleanup_function_t));

    flint_cleanup_functions[flint_num_cleanup_functions] = cleanup_function;

    flint_num_cleanup_functions++;

#if FLINT_REENTRANT && !HAVE_TLS
    pthread_mutex_unlock(&register_lock);
#endif
}

void _fmpz_cleanup();

void flint_cleanup()
{
    size_t i;

#if FLINT_REENTRANT && !HAVE_TLS
    pthread_mutex_lock(&register_lock);
#endif

    for (i = 0; i < flint_num_cleanup_functions; i++)
        flint_cleanup_functions[i]();

    flint_free(flint_cleanup_functions);
    flint_cleanup_functions = NULL;
    flint_num_cleanup_functions = 0;

    mpfr_free_cache();
    _fmpz_cleanup();

#if FLINT_REENTRANT && !HAVE_TLS
    pthread_mutex_unlock(&register_lock);
#endif

}


