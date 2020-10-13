/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"

void
_ca_vec_clear(ca_ptr v, slong n, ca_ctx_t ctx)
{
    slong i;
    for (i = 0; i < n; i++)
        ca_clear(v + i, ctx);
    flint_free(v);
}

void
ca_vec_clear(ca_vec_t vec, ca_ctx_t ctx)
{
    if (vec->entries != NULL)
    {
        _ca_vec_clear(vec->entries, vec->alloc, ctx);
    }
}
