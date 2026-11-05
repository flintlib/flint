/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "impl.h"

typedef struct
{
    int type;
    void * obj;
    _acb_mzeta_cache_clear_func clear;
}
cache_entry;

void
acb_mzeta_cache_init(acb_mzeta_cache_t cache)
{
    cache->entries = NULL;
    cache->length = 0;
    cache->alloc = 0;
}

void
acb_mzeta_cache_clear(acb_mzeta_cache_t cache)
{
    cache_entry * e = cache->entries;
    slong i;

    for (i = 0; i < cache->length; i++)
    {
        e[i].clear(e[i].obj);
        flint_free(e[i].obj);
    }
    flint_free(cache->entries);
}

void *
_acb_mzeta_cache_find(acb_mzeta_cache_t cache, int type,
    _acb_mzeta_cache_match_func match, const void * key)
{
    cache_entry * e = cache->entries;
    slong i;

    /* most recent first */
    for (i = cache->length - 1; i >= 0; i--)
        if (e[i].type == type && match(e[i].obj, key))
            return e[i].obj;

    return NULL;
}

/* takes ownership of obj (allocated with flint_malloc) */
void
_acb_mzeta_cache_insert(acb_mzeta_cache_t cache, int type, void * obj,
    _acb_mzeta_cache_clear_func clear)
{
    cache_entry * e = cache->entries;
    slong i, count = 0, oldest = -1;

    for (i = 0; i < cache->length; i++)
    {
        if (e[i].type == type)
        {
            if (oldest == -1)
                oldest = i;
            count++;
        }
    }

    if (count >= _ACB_MZETA_CACHE_MAX_PER_TYPE)
    {
        e[oldest].clear(e[oldest].obj);
        flint_free(e[oldest].obj);
        for (i = oldest; i < cache->length - 1; i++)
            e[i] = e[i + 1];
        cache->length--;
    }

    if (cache->length == cache->alloc)
    {
        cache->alloc = FLINT_MAX(4, 2 * cache->alloc);
        cache->entries = flint_realloc(cache->entries, cache->alloc * sizeof(cache_entry));
        e = cache->entries;
    }

    e[cache->length].type = type;
    e[cache->length].obj = obj;
    e[cache->length].clear = clear;
    cache->length++;
}
