/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"

/* The current version is published with a release store and read with
   an acquire load, so that a reader seeing a version also sees its
   contents; everything else happens under the cache's lock.  Without
   pthreads there is one thread and plain accesses do. */
#if FLINT_USES_PTHREAD
# if defined(__GNUC__)
#  define LOAD_ACQUIRE(p) __atomic_load_n(p, __ATOMIC_ACQUIRE)
#  define STORE_RELEASE(p, v) __atomic_store_n(p, v, __ATOMIC_RELEASE)
# else
#  include <stdatomic.h>
/* the field is only ever accessed through these two macros */
#  define LOAD_ACQUIRE(p) \
    atomic_load_explicit((_Atomic(void *) *) (p), memory_order_acquire)
#  define STORE_RELEASE(p, v) \
    atomic_store_explicit((_Atomic(void *) *) (p), v, memory_order_release)
# endif
#else
# define LOAD_ACQUIRE(p) (*(p))
# define STORE_RELEASE(p, v) (*(p) = (v))
#endif

/* a retired version */
typedef struct _retired_struct
{
    void * version;
    struct _retired_struct * next;
}
_retired_struct;

/* the registry of caches that have published a version, for
   flint_cleanup_master */
static flint_shared_cache_struct * _registry = NULL;
#if FLINT_USES_PTHREAD
static pthread_mutex_t _registry_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

void *
flint_shared_cache_get(flint_shared_cache_struct * cache)
{
    return LOAD_ACQUIRE(&cache->cur);
}

void
flint_shared_cache_lock(flint_shared_cache_struct * cache)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&cache->lock);
#else
    (void) cache;
#endif
}

void
flint_shared_cache_unlock(flint_shared_cache_struct * cache)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&cache->lock);
#else
    (void) cache;
#endif
}

/* under the cache's lock: version replaces the current one, which is
   retired */
void
flint_shared_cache_publish(flint_shared_cache_struct * cache, void * version)
{
    void * old = LOAD_ACQUIRE(&cache->cur);

    if (old != NULL)
    {
        _retired_struct * r = flint_malloc(sizeof(_retired_struct));
        r->version = old;
        r->next = cache->retired;
        cache->retired = r;
    }

    STORE_RELEASE(&cache->cur, version);

    if (!cache->registered)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&_registry_lock);
#endif
        cache->next = _registry;
        _registry = cache;
        cache->registered = 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&_registry_lock);
#endif
    }
}

/* the current version if sufficient(version, args), else the version
   build(current, args), published; NULL is never sufficient */
void *
flint_shared_cache_ensure(flint_shared_cache_struct * cache,
    flint_shared_cache_sufficient_t sufficient,
    flint_shared_cache_build_t build, void * args)
{
    void * v = LOAD_ACQUIRE(&cache->cur);

    if (v != NULL && sufficient(v, args))
        return v;

    flint_shared_cache_lock(cache);
    v = LOAD_ACQUIRE(&cache->cur);
    if (v == NULL || !sufficient(v, args))
    {
        v = build(v, args);
        flint_shared_cache_publish(cache, v);
    }
    flint_shared_cache_unlock(cache);

    return v;
}

/* under the cache's lock, or with the registry's: free every version */
static void
_flint_shared_cache_free_all(flint_shared_cache_struct * cache)
{
    void * v = LOAD_ACQUIRE(&cache->cur);
    _retired_struct * r, * next;

    STORE_RELEASE(&cache->cur, NULL);
    if (v != NULL)
        cache->free_version(v);

    for (r = cache->retired; r != NULL; r = next)
    {
        next = r->next;
        cache->free_version(r->version);
        flint_free(r);
    }
    cache->retired = NULL;
}

/* frees every version of the cache: only when no other thread can be
   reading it (tuning programs measuring cold rebuilds) */
void
flint_shared_cache_clear(flint_shared_cache_struct * cache)
{
    flint_shared_cache_lock(cache);
    _flint_shared_cache_free_all(cache);
    flint_shared_cache_unlock(cache);
}

/* the number of caches that have published a version since the last
   flint_cleanup_master (tests checking that a regime touches no cache) */
slong
_flint_shared_cache_count(void)
{
    flint_shared_cache_struct * cache;
    slong count = 0;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&_registry_lock);
#endif
    for (cache = _registry; cache != NULL; cache = cache->next)
        count++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&_registry_lock);
#endif
    return count;
}

/* flint_cleanup_master: frees every version of every cache that has
   published one */
void
_flint_shared_cache_cleanup_all(void)
{
    flint_shared_cache_struct * list, * cache, * next;

    /* the registry is detached first: publish takes the locks in the
       other order */
#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&_registry_lock);
#endif
    list = _registry;
    _registry = NULL;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&_registry_lock);
#endif

    for (cache = list; cache != NULL; cache = next)
    {
        flint_shared_cache_lock(cache);
        next = cache->next;
        _flint_shared_cache_free_all(cache);
        cache->next = NULL;
        cache->registered = 0;
        flint_shared_cache_unlock(cache);
    }
}
