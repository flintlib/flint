/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"

ca_ptr
_ca_vec_init(slong n, ca_ctx_t ctx)
{
    slong i;
    ca_ptr v = (ca_ptr) flint_malloc(sizeof(ca_struct) * n);

    for (i = 0; i < n; i++)
        ca_init(v + i, ctx);

    return v;
}

void
ca_vec_init(ca_vec_t res, slong len, ca_ctx_t ctx)
{
    if (len == 0)
    {
        res->entries = NULL;
        res->length = 0;
        res->alloc = 0;
    }
    else
    {
        res->entries = _ca_vec_init(len, ctx);
        res->length = len;
        res->alloc = len;
    }
}
