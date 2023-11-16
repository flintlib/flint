/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"

#define INITIAL_HASH_SIZE 16

void
ca_ext_cache_init(ca_ext_cache_t cache, ca_ctx_t ctx)
{
    slong i;

    cache->items = NULL;
    cache->length = 0;
    cache->alloc = 0;

    cache->hash_size = INITIAL_HASH_SIZE;
    cache->hash_table = flint_malloc(sizeof(slong) * INITIAL_HASH_SIZE);
    for (i = 0; i < INITIAL_HASH_SIZE; i++)
        cache->hash_table[i] = -1;
}

