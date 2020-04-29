/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "flint.h"
#include "hashmap.h"

void hashmap_init(hashmap_t h, slong size)
{
    memset(h, 0, sizeof(*h));
    h->size = (UWORD(1)) << MIN_HASHMAP_BITS;
    while (h->size < size) h->size <<= UWORD(1);

    h->mask = 2*h->size - 1;
    h->num = 0;
    h->table = flint_calloc(2*h->size, sizeof(*h->table));
    h->keys = flint_malloc(h->size*sizeof(*h->keys));
    h->vals = flint_malloc(h->size*sizeof(*h->vals));
}

void hashmap_clear(hashmap_t h)
{
    flint_free(h->table);
    flint_free(h->keys);
    flint_free(h->vals);
    memset(h, 0, sizeof(*h));
}

static void _hashmap_rehash(hashmap_t h)
{
    slong i, num = h->num;
    h->size <<= 1;
    h->mask = 2*h->size - 1;
    h->table = realloc(h->table, 2*h->size*sizeof(*h->table));
    h->keys = realloc(h->keys, h->size*sizeof(*h->keys));
    h->vals = realloc(h->vals, h->size*sizeof(*h->vals));
    memset(h->table, 0, 2*h->size*sizeof(*h->table));
    h->num = 0;
    for (i = 0; i < num; ++i)
        hashmap_put(h, h->keys[i], h->vals[i]);
}

static slong _hashmap_pos(hashmap_t h, slong key, int skip_deleted)
{
    slong ind, pos = key*UWORD(13282407956253574709) + UWORD(286824421);
    for(; (ind = h->table[pos & h->mask]); pos++)
    {
        if (!skip_deleted && ind == -1) break;
        if (ind != -1 && h->keys[ind - 1] == key) break;
    }
    return pos & h->mask;
}

void * hashmap_get(hashmap_t h, slong key)
{
    slong pos = _hashmap_pos(h, key, 1), ind = h->table[pos];
    return (ind = h->table[pos]) ? h->vals[ind - 1] : NULL;
}

void hashmap_put(hashmap_t h, slong key, void *val)
{
    slong pos = _hashmap_pos(h, key, 1), ind = h->table[pos];
    if(ind > 0) h->vals[ind - 1] = val;
    else if(h->num < h->size)
    {
        pos = _hashmap_pos(h, key, 0);
        h->keys[h->num] = key;
        h->vals[h->num] = val;
        h->table[pos] = ++h->num;
    }
    else
    {
        _hashmap_rehash(h);
        hashmap_put(h, key, val);
    }
}

void hashmap_rem(hashmap_t h, slong key)
{
    slong pos = _hashmap_pos(h, key, 1), ind = h->table[pos];
    if(ind > 0)
    {
        h->table[pos] = -1; /* Mark as deleted, new entries can be put here */
        if (ind < h->num)
        {
            h->keys[ind - 1] = h->keys[h->num - 1];
            h->vals[ind - 1] = h->vals[h->num - 1];
            h->table[_hashmap_pos(h, h->keys[ind - 1], 1)] = ind;
        }
        h->num--;
    }
}