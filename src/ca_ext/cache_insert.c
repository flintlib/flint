/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"

ca_ext_ptr ca_ext_cache_insert(ca_ext_cache_t cache, const ca_ext_t x, ca_ctx_t ctx)
{
    ulong xhash;
    slong i, loc;

    xhash = ca_ext_hash(x, ctx);

    /* make room for inserting entry if needed */
    if (cache->length == cache->alloc)
    {
        slong new_alloc;

        new_alloc = FLINT_MAX(1, cache->alloc * 2);
        cache->items = flint_realloc(cache->items, sizeof(ca_ext_struct *) * new_alloc);

        for (i = cache->alloc; i < new_alloc; i++)
            cache->items[i] = flint_malloc(sizeof(ca_ext_struct));

        cache->alloc = new_alloc;
    }

    /* rehash if needed */
    if (cache->length >= 0.5 * cache->hash_size)
    {
        slong new_size, j;
        slong * new_table;
        ulong thash;

        new_size = cache->hash_size * 2;

        new_table = flint_malloc(sizeof(slong) * new_size);

        for (i = 0; i < new_size; i++)
            new_table[i] = -1;

        for (i = 0; i < cache->length; i++)
        {
            thash = ca_ext_hash(cache->items[i], ctx);

            j = thash % new_size;

            while (new_table[j] != -1)
            {
                j++;
                if (j == new_size)
                    j = 0;
            }

            new_table[j] = i;
        }

        flint_free(cache->hash_table);

        cache->hash_size = new_size;
        cache->hash_table = new_table;
    }

    loc = xhash % ((ulong) cache->hash_size);

    for (i = 0; i < cache->hash_size; i++)
    {
        /* not found, so insert */
        if (cache->hash_table[loc] == -1)
        {
            ca_ext_init_set(cache->items[cache->length], x, ctx);
            cache->hash_table[loc] = cache->length;
            cache->length++;
            return cache->items[cache->length - 1];
        }

        /* found */
        if (ca_ext_equal_repr(cache->items[cache->hash_table[loc]], x, ctx))
            return cache->items[cache->hash_table[loc]];

        loc++;
        if (loc == cache->hash_size)
            loc = 0;
    }

    /* cannot happen */
    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
}

