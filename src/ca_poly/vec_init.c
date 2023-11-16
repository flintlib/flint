/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

ca_poly_struct *
_ca_poly_vec_init(slong n, ca_ctx_t ctx)
{
    slong i;
    ca_poly_struct * v = (ca_poly_struct *) flint_malloc(sizeof(ca_poly_struct) * n);

    for (i = 0; i < n; i++)
        ca_poly_init(v + i, ctx);

    return v;
}

void
ca_poly_vec_init(ca_poly_vec_t res, slong len, ca_ctx_t ctx)
{
    if (len == 0)
    {
        res->entries = NULL;
        res->length = 0;
        res->alloc = 0;
    }
    else
    {
        res->entries = _ca_poly_vec_init(len, ctx);
        res->length = len;
        res->alloc = len;
    }
}
