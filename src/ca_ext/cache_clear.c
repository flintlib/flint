/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"

void
ca_ext_cache_clear(ca_ext_cache_t cache, ca_ctx_t ctx)
{
    slong i;

    /* reverse order ensures that we free f(x) before freeing any element
       used by x */
    for (i = cache->length - 1; i >= 0; i--)
        ca_ext_clear(cache->items[i], ctx);

    for (i = 0; i < cache->alloc; i++)
        flint_free(cache->items[i]);

    flint_free(cache->items);
    flint_free(cache->hash_table);
}
