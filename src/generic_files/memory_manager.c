/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2016 Claus Fieker
    Copyright (C) 2016 William Hart.
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "mpfr.h"
#include "thread_pool.h"

#if FLINT_USES_GC
#include "gc.h"
#endif

static void * _flint_malloc(size_t);
static void * _flint_calloc(size_t, size_t);
static void * _flint_realloc(void *, size_t);
static void _flint_free(void *);

static void  *(*__flint_allocate_func) (size_t) = _flint_malloc;
static void  *(*__flint_callocate_func) (size_t, size_t) = _flint_calloc;
static void  *(*__flint_reallocate_func) (void *, size_t) = _flint_realloc;
static void  (*__flint_free_func) (void *) = _flint_free;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
#include <pthread.h>

static pthread_once_t register_initialised = PTHREAD_ONCE_INIT;
pthread_mutex_t register_lock;

static pthread_once_t alloc_func_init = PTHREAD_ONCE_INIT;
pthread_mutex_t alloc_func_lock;

void __flint_set_memory_functions_init(void)
{
   pthread_mutex_init(&alloc_func_lock, NULL);
}
#endif

static void flint_memory_error(size_t size)
{
    flint_throw(FLINT_ERROR, "Unable to allocate memory (%zu).\n", size);
}

void __flint_get_memory_functions(void *(**alloc_func) (size_t),
                             void *(**calloc_func) (size_t, size_t),
                             void *(**realloc_func) (void *, size_t),
                             void (**free_func) (void *))
{
#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&alloc_func_init, __flint_set_memory_functions_init);
    pthread_mutex_lock(&alloc_func_lock);
#endif

    *alloc_func = __flint_allocate_func;
    *calloc_func = __flint_callocate_func;
    *realloc_func = __flint_reallocate_func;
    *free_func = __flint_free_func;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&alloc_func_lock);
#endif
}

void __flint_set_memory_functions(void *(*alloc_func) (size_t),
                             void *(*calloc_func) (size_t, size_t),
                             void *(*realloc_func) (void *, size_t),
                             void (*free_func) (void *))
{
#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&alloc_func_init, __flint_set_memory_functions_init);
    pthread_mutex_lock(&alloc_func_lock);
#endif

  __flint_allocate_func = alloc_func;
  __flint_callocate_func = calloc_func;
  __flint_reallocate_func = realloc_func;
  __flint_free_func = free_func;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&alloc_func_lock);
#endif
}

void * _flint_malloc(size_t size)
{
   void * ptr;

#if FLINT_USES_GC
   ptr = GC_malloc(size);
#else
   ptr = malloc(size);
#endif

   return ptr;
}

FLINT_WARN_UNUSED void * flint_malloc(size_t size)
{
   void * ptr = (*__flint_allocate_func)(size);

   if (ptr == NULL)
        flint_memory_error(size);

   return ptr;
}

void * _flint_realloc(void * ptr, size_t size)
{
    void * ptr2;

#if FLINT_USES_GC
    ptr2 = GC_realloc(ptr, size);
#else
    ptr2 = realloc(ptr, size);
#endif

    return ptr2;
}


FLINT_WARN_UNUSED void * flint_realloc(void * ptr, size_t size)
{
    void * ptr2;

    if (ptr)
      ptr2 = (*__flint_reallocate_func)(ptr, size);
    else
      ptr2 = (*__flint_allocate_func)(size);

    if (ptr2 == NULL)
        flint_memory_error(size);

    return ptr2;
}

void * _flint_calloc(size_t num, size_t size)
{
   void * ptr;

#if FLINT_USES_GC
    ptr = GC_malloc(num*size);
#else
    ptr = calloc(num, size);
#endif

    return ptr;
}

FLINT_WARN_UNUSED void * flint_calloc(size_t num, size_t size)
{
   void * ptr;

    ptr = (*__flint_callocate_func)(num, size);

    if (ptr == NULL)
        flint_memory_error(size);

    return ptr;
}

void _flint_free(void * ptr)
{
#if !FLINT_USES_GC
    free(ptr);
#endif
}

void flint_free(void * ptr)
{
   (*__flint_free_func)(ptr);
}


FLINT_TLS_PREFIX size_t flint_num_cleanup_functions = 0;

FLINT_TLS_PREFIX flint_cleanup_function_t * flint_cleanup_functions = NULL;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
void register_init(void)
{
   pthread_mutex_init(&register_lock, NULL);
}
#endif

void flint_register_cleanup_function(flint_cleanup_function_t cleanup_function)
{
#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&register_initialised, register_init);
    pthread_mutex_lock(&register_lock);
#endif

    flint_cleanup_functions = flint_realloc(flint_cleanup_functions,
        (flint_num_cleanup_functions + 1) * sizeof(flint_cleanup_function_t));

    flint_cleanup_functions[flint_num_cleanup_functions] = cleanup_function;

    flint_num_cleanup_functions++;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&register_lock);
#endif
}

void _fmpz_cleanup(void);

void _flint_cleanup(void)
{
    size_t i;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_lock(&register_lock);
#endif

    for (i = 0; i < flint_num_cleanup_functions; i++)
        flint_cleanup_functions[i]();

    flint_free(flint_cleanup_functions);
    flint_cleanup_functions = NULL;
    flint_num_cleanup_functions = 0;

    mpfr_free_cache();
    _fmpz_cleanup();

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&register_lock);
#endif

}

void flint_cleanup(void)
{
#if !FLINT_REENTRANT || FLINT_USES_TLS
   _flint_cleanup();
#endif
}

void flint_cleanup_master(void)
{
    if (global_thread_pool_initialized)
    {
        thread_pool_clear(global_thread_pool);
        global_thread_pool_initialized = 0;
    }
    _flint_cleanup();
}
