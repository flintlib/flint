/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_field.h"

void
ca_field_cache_clear(ca_field_cache_t cache, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < cache->length; i++)
        ca_field_clear(cache->items[i], ctx);

    for (i = 0; i < cache->alloc; i++)
        flint_free(cache->items[i]);

    flint_free(cache->items);
    flint_free(cache->hash_table);
}
