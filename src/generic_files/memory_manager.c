/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2016 Claus Fieker
    Copyright (C) 2016 William Hart.
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <mpfr.h>
#include "thread_pool.h"

#if FLINT_USES_GC
# include <gc.h>
#endif

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
# include <pthread.h>
#endif

#if defined(_MSC_VER) || defined(__MINGW32__) || defined(__MINGW64__)
# include <malloc.h>
#endif

/* memory functions **********************************************************/

static void * _flint_malloc(size_t);
static void * _flint_calloc(size_t, size_t);
static void * _flint_realloc(void *, size_t);
static void _flint_free(void *);
static void * _flint_aligned_alloc(size_t, size_t);
static void * _flint_aligned_alloc2(size_t, size_t);
static void _flint_aligned_free(void *);
static void _flint_aligned_free2(void *);

static void * (* __flint_allocate_func)(size_t) = _flint_malloc;
static void * (* __flint_callocate_func)(size_t, size_t) = _flint_calloc;
static void * (* __flint_reallocate_func)(void *, size_t) = _flint_realloc;
static void (* __flint_free_func)(void *) = _flint_free;
#if HAVE_ALIGNED_ALLOC || HAVE__ALIGNED_MALLOC
static void * (* __flint_aligned_allocate_func)(size_t, size_t) = _flint_aligned_alloc;
static void (* __flint_aligned_free_func)(void *) = _flint_aligned_free;
#else
static void * (* __flint_aligned_allocate_func)(size_t, size_t) = _flint_aligned_alloc2;
static void (* __flint_aligned_free_func)(void *) = _flint_aligned_free2;
#endif

FLINT_STATIC_NOINLINE void flint_memory_error(size_t size)
{
    flint_throw(FLINT_ERROR, "Unable to allocate memory (%zu).\n", size);
}

FLINT_STATIC_NOINLINE void flint_aligned_memory_error(size_t alignment, size_t size)
{
    flint_throw(FLINT_ERROR, "Unable to allocate %zu bytes with alignment %zu\n", size, alignment);
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

void * _flint_aligned_alloc(size_t alignment, size_t size)
{
#if HAVE__ALIGNED_MALLOC
    return _aligned_malloc(size, alignment);
#elif HAVE_ALIGNED_ALLOC
    return aligned_alloc(alignment, size);
#else
    return NULL;
#endif
}

/* NOTE: This function assumes alignment >= sizeof(ulong) */
void * _flint_aligned_alloc2(size_t alignment, size_t size)
{
    size_t alloc_size;
    void * alloc_ptr;
    ulong * ret_ptr;

    alloc_size = size + alignment;

    alloc_ptr = flint_malloc(alloc_size);

    /* Case 1: alloc_ptr aligned with (alignment, alignment - sizeof(ulong)).
               We only need `size + sizeof(ulong)' bytes.

       Case 2: alloc_ptr aligned with (alignment, n),
               where 0 <= n < alignment - sizeof(ulong). We will not use the
               first `alignment - n - sizeof(ulong)' bytes, and hence we
               need `size + alignment - n - sizeof(ulong)' bytes. */

    ret_ptr = (ulong *) ((((size_t) alloc_ptr) & (-alignment)) + alignment);
    ret_ptr[-1] = (char *) ret_ptr - (char *) alloc_ptr;

    return ret_ptr;
}

FLINT_WARN_UNUSED void * flint_aligned_alloc(size_t alignment, size_t size)
{
    void * ptr;

    FLINT_ASSERT(size % alignment == 0);

    ptr = (*__flint_aligned_allocate_func)(alignment, size);

    if (ptr == NULL)
        flint_aligned_memory_error(alignment, size);

    return ptr;
}

void _flint_aligned_free(void * p)
{
#if HAVE__ALIGNED_MALLOC
    _aligned_free(p);
#elif HAVE_ALIGNED_ALLOC
    free(p);
#else
    return;
#endif
}

void _flint_aligned_free2(void * p)
{
    size_t * ptr = p;
    if (ptr != NULL)
        flint_free((char *) ptr - ptr[-1]);
}

void flint_aligned_free(void * ptr)
{
    (*__flint_aligned_free_func)(ptr);
}

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
static pthread_once_t register_initialised = PTHREAD_ONCE_INIT;
pthread_mutex_t register_lock;

static pthread_once_t alloc_func_init = PTHREAD_ONCE_INIT;
pthread_mutex_t alloc_func_lock;

void __flint_set_memory_functions_init(void)
{
    pthread_mutex_init(&alloc_func_lock, NULL);
}
#endif

/* aligned memory getter and setter ******************************************/

void __flint_get_all_memory_functions(
        void * (** alloc_func)(size_t),
        void * (** calloc_func)(size_t, size_t),
        void * (** realloc_func)(void *, size_t),
        void (** free_func)(void *),
        void * (** aligned_alloc_func)(size_t, size_t),
        void (** aligned_free_func)(void *))
{
#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&alloc_func_init, __flint_set_memory_functions_init);
    pthread_mutex_lock(&alloc_func_lock);
#endif

    *alloc_func = __flint_allocate_func;
    *calloc_func = __flint_callocate_func;
    *realloc_func = __flint_reallocate_func;
    *free_func = __flint_free_func;
    *aligned_alloc_func = __flint_aligned_allocate_func;
    *aligned_free_func = __flint_aligned_free_func;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&alloc_func_lock);
#endif
}

void __flint_set_all_memory_functions(
        void * (* alloc_func)(size_t),
        void * (* calloc_func)(size_t, size_t),
        void * (* realloc_func)(void *, size_t),
        void (* free_func)(void *),
        void * (* aligned_alloc_func)(size_t, size_t),
        void (* aligned_free_func)(void *))
{
#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&alloc_func_init, __flint_set_memory_functions_init);
    pthread_mutex_lock(&alloc_func_lock);
#endif

    __flint_allocate_func = alloc_func;
    __flint_callocate_func = calloc_func;
    __flint_reallocate_func = realloc_func;
    __flint_free_func = free_func;
    __flint_aligned_allocate_func = aligned_alloc_func;
    __flint_aligned_free_func = aligned_free_func;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&alloc_func_lock);
#endif
}

/* non-aligned memory getter and setter **************************************/

void __flint_get_memory_functions(
        void * (** alloc_func)(size_t),
        void * (** calloc_func)(size_t, size_t),
        void * (** realloc_func)(void *, size_t),
        void (** free_func)(void *))
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

void __flint_set_memory_functions(
        void * (* alloc_func)(size_t),
        void * (* calloc_func)(size_t, size_t),
        void * (* realloc_func)(void *, size_t),
        void (* free_func)(void *))
{
#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_once(&alloc_func_init, __flint_set_memory_functions_init);
    pthread_mutex_lock(&alloc_func_lock);
#endif

    __flint_allocate_func = alloc_func;
    __flint_callocate_func = calloc_func;
    __flint_reallocate_func = realloc_func;
    __flint_free_func = free_func;

    __flint_aligned_allocate_func = _flint_aligned_alloc2;
    __flint_aligned_free_func = _flint_aligned_free2;

#if FLINT_REENTRANT && !FLINT_USES_TLS && FLINT_USES_PTHREAD
    pthread_mutex_unlock(&alloc_func_lock);
#endif
}

/* cleanup functions *********************************************************/

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

static void _flint_cleanup(void)
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
