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
_ca_poly_vec_fit_length(ca_poly_vec_t vec, slong len, ca_ctx_t ctx)
{
    if (len > vec->alloc)
    {
        slong i;

        if (len < 2 * vec->alloc)
            len = 2 * vec->alloc;

        vec->entries = flint_realloc(vec->entries, len * sizeof(ca_poly_struct));

        for (i = vec->alloc; i < len; i++)
            ca_poly_init(vec->entries + i, ctx);

        vec->alloc = len;
    }
}

void
ca_poly_vec_set_length(ca_poly_vec_t vec, slong len, ca_ctx_t ctx)
{
    slong i;

    if (vec->length > len)
    {
        for (i = len; i < vec->length; i++)
            ca_poly_zero(vec->entries + i, ctx);
    }
    else if (vec->length < len)
    {
        _ca_poly_vec_fit_length(vec, len, ctx);

        for (i = vec->length; i < len; i++)
            ca_poly_zero(vec->entries + i, ctx);
    }

    vec->length = len;
}
