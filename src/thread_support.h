/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef THREAD_SUPPORT_H
#define THREAD_SUPPORT_H

#include "flint.h"

#if FLINT_USES_PTHREAD
# include <pthread.h>
#endif

#ifdef __cplusplus
 extern "C" {
#endif

#define FLINT_DEFAULT_THREAD_LIMIT 99999

slong flint_request_threads(thread_pool_handle ** handles, slong thread_limit);

void flint_give_back_threads(thread_pool_handle * handles, slong num_handles);

slong flint_get_num_available_threads(void);

#define FLINT_PARALLEL_UNIFORM 1
#define FLINT_PARALLEL_STRIDED 2
#define FLINT_PARALLEL_DYNAMIC 4
#define FLINT_PARALLEL_BSPLIT_LEFT_INPLACE 8
#define FLINT_PARALLEL_VERBOSE 512

typedef void (* do_func_t)(slong i, void * args);

void flint_parallel_do(do_func_t f, void * args, slong n, int thread_limit, int flags);

typedef void (* bsplit_merge_func_t)(void *, void *, void *, void *);
typedef void (* bsplit_basecase_func_t)(void *, slong, slong, void *);
typedef void (* bsplit_init_func_t)(void *, void *);
typedef void (* bsplit_clear_func_t)(void *, void *);

void flint_parallel_binary_splitting(void * res, bsplit_basecase_func_t basecase, bsplit_merge_func_t merge,
    size_t sizeof_res, bsplit_init_func_t init, bsplit_clear_func_t clear, void * args, slong a, slong b, slong basecase_cutoff, int thread_limit, int flags);

/* Shared caches **************************************************************/

/* A shared cache holds one process-wide, read-only "version" of some
   lazily computed data -- a table of constants, say -- that any thread
   reads without locking (flint_shared_cache_get: one acquire load).
   It grows by replacement: a thread finding the current version
   insufficient takes the lock, checks again (another thread may have
   grown the cache meanwhile), builds a larger version, usually from the
   current one, and publishes it.  Superseded versions are retired, not
   freed, as other threads may still be reading them, so the versions of
   a cache should grow geometrically; every version is freed by
   flint_cleanup_master() (NOT by the per-thread flint_cleanup()).  The
   whole protocol is flint_shared_cache_ensure. */

typedef void (* flint_shared_cache_free_t)(void * version);

typedef struct flint_shared_cache_struct
{
#if FLINT_USES_PTHREAD
    pthread_mutex_t lock;
#endif
    void * cur;                     /* the current version, or NULL
                                       (accessed atomically) */
    void * retired;                 /* superseded versions */
    flint_shared_cache_free_t free_version;
    struct flint_shared_cache_struct * next;    /* the cleanup registry */
    int registered;
}
flint_shared_cache_struct;

/* static initializer: free_version frees one version */
#if FLINT_USES_PTHREAD
# define FLINT_SHARED_CACHE_INIT(free_version) \
    { PTHREAD_MUTEX_INITIALIZER, NULL, NULL, free_version, NULL, 0 }
#else
# define FLINT_SHARED_CACHE_INIT(free_version) \
    { NULL, NULL, free_version, NULL, 0 }
#endif

void * flint_shared_cache_get(flint_shared_cache_struct * cache);

void flint_shared_cache_lock(flint_shared_cache_struct * cache);
void flint_shared_cache_unlock(flint_shared_cache_struct * cache);

void flint_shared_cache_publish(flint_shared_cache_struct * cache,
    void * version);

typedef int (* flint_shared_cache_sufficient_t)(const void * version,
    void * args);
typedef void * (* flint_shared_cache_build_t)(const void * version,
    void * args);

void * flint_shared_cache_ensure(flint_shared_cache_struct * cache,
    flint_shared_cache_sufficient_t sufficient,
    flint_shared_cache_build_t build, void * args);

void flint_shared_cache_clear(flint_shared_cache_struct * cache);

void _flint_shared_cache_cleanup_all(void);
slong _flint_shared_cache_count(void);

#ifdef __cplusplus
}
#endif

#endif
