/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"
#include "ca_field.h"

ca_field_ptr ca_field_cache_lookup_qqbar(ca_field_cache_t cache, const qqbar_t x, ca_ctx_t ctx)
{
    ulong xhash;
    ca_field_ptr K;
    slong i, loc;

    xhash = qqbar_hash(x);

    loc = xhash % ((ulong) cache->hash_size);

    for (i = 0; i < cache->hash_size; i++)
    {
        /* not found */
        if (cache->hash_table[loc] == -1)
            return NULL;

        K = cache->items[cache->hash_table[loc]];
        /* found */
        if (CA_FIELD_IS_NF(K) && qqbar_equal(x, CA_FIELD_NF_QQBAR(K)))
            return K;

        loc++;
        if (loc == cache->hash_size)
            loc = 0;
    }

    /* cannot happen */
    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
}

ulong
_ca_field_hash(ca_ext_struct ** ext, slong len, ca_ctx_t ctx)
{
    ulong s;
    slong i;

    s = 0;
    for (i = 0; i < len; i++)
        s = s * CA_FIELD_HASH_C + CA_EXT_HASH(ext[i]);

    return s;
}

ulong
ca_field_hash(const ca_field_t K, ca_ctx_t ctx)
{
    return _ca_field_hash(K->ext, K->length, ctx);
}

static int _ca_field_equal_ext(const ca_field_t K, ca_ext_struct ** x, slong len, ca_ctx_t ctx)
{
    slong i;

    if (len != CA_FIELD_LENGTH(K))
        return 0;

    for (i = 0; i < len; i++)
        if (CA_FIELD_EXT_ELEM(K, i) != x[i])
            return 0;

    return 1;
}

void
ca_field_init_set_ext(ca_field_t K, ca_ext_struct ** ext, slong len, ca_ctx_t ctx)
{
    if (len == 0)
    {
        ca_field_init_qq(K, ctx);
    }
    else if (len == 1 && CA_EXT_IS_QQBAR(ext[0]))
    {
        CA_FIELD_LENGTH(K) = 1;
        CA_FIELD_EXT(K) = flint_malloc(sizeof(ca_ext_ptr));
        CA_FIELD_EXT_ELEM(K, 0) = ext[0];
        CA_FIELD_IDEAL_P(K) = NULL;
        CA_FIELD_IDEAL_LENGTH(K) = -1;
        CA_FIELD_IDEAL_ALLOC(K) = 0;
        CA_FIELD_HASH(K) = CA_EXT_HASH(ext[0]);
    }
    else
    {
        slong i;

        ca_field_init_multi(K, len, ctx);

        for (i = 0; i < len; i++)
            ca_field_set_ext(K, i, ext[i], ctx);
    }
}


ca_field_ptr ca_field_cache_insert_ext(ca_field_cache_t cache, ca_ext_struct ** x, slong length, ca_ctx_t ctx)
{
    ulong xhash;
    slong i, loc;

    xhash = _ca_field_hash(x, length, ctx);

    /* make room for inserting entry if needed */
    if (cache->length == cache->alloc)
    {
        slong new_alloc;

        new_alloc = FLINT_MAX(1, cache->alloc * 2);
        cache->items = flint_realloc(cache->items, sizeof(ca_field_struct *) * new_alloc);

        for (i = cache->alloc; i < new_alloc; i++)
            cache->items[i] = flint_malloc(sizeof(ca_field_struct));

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
            thash = ca_field_hash(cache->items[i], ctx);

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
            ca_field_ptr res;

            ca_field_init_set_ext(cache->items[cache->length], x, length, ctx);
            cache->hash_table[loc] = cache->length;
            cache->length++;

            /* save pointer; build_ideal can resize the cache */
            res = cache->items[cache->length - 1];

            ca_field_build_ideal(res, ctx);

            return res;
        }

        /* found */
        if (_ca_field_equal_ext(cache->items[cache->hash_table[loc]], x, length, ctx))
            return cache->items[cache->hash_table[loc]];

        loc++;
        if (loc == cache->hash_size)
            loc = 0;
    }

    /* cannot happen */
    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
}
