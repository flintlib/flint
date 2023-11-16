/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_vec_clear(ca_poly_struct * v, slong n, ca_ctx_t ctx)
{
    slong i;
    for (i = 0; i < n; i++)
        ca_poly_clear(v + i, ctx);
    flint_free(v);
}

void
ca_poly_vec_clear(ca_poly_vec_t vec, ca_ctx_t ctx)
{
    if (vec->entries != NULL)
    {
        _ca_poly_vec_clear(vec->entries, vec->alloc, ctx);
    }
}
