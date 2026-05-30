/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

void
gr_poly_vec_init(gr_poly_vec_t vec, slong len, gr_ctx_t ctx)
{
    slong i;

    vec->length = vec->alloc = len;

    if (len == 0)
    {
        vec->entries = NULL;
    }
    else
    {
        vec->entries = flint_malloc(len * sizeof(gr_poly_struct));
        for (i = 0; i < len; i++)
            gr_poly_init(vec->entries + i, ctx);
    }
}

void
gr_poly_vec_clear(gr_poly_vec_t vec, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < vec->alloc; i++)
        gr_poly_clear(vec->entries + i, ctx);
    flint_free(vec->entries);
}

int
gr_poly_vec_set(gr_poly_vec_t res, const gr_poly_vec_t src, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;

    if (res == src)
        return GR_SUCCESS;

    gr_poly_vec_set_length(res, src->length, ctx);
    for (i = 0; i < src->length; i++)
        status |= gr_poly_set(res->entries + i, src->entries + i, ctx);
    return status;
}

void
gr_poly_vec_fit_length(gr_poly_vec_t vec, slong len, gr_ctx_t ctx)
{
    slong i, alloc = vec->alloc;

    if (len > alloc)
    {
        if (len < 2 * alloc)
            len = 2 * alloc;

        vec->entries = flint_realloc(vec->entries, len * sizeof(gr_poly_struct));
        for (i = alloc; i < len; i++)
            gr_poly_init(vec->entries + i, ctx);
        vec->alloc = len;
    }
}

void
gr_poly_vec_set_length(gr_poly_vec_t vec, slong len, gr_ctx_t ctx)
{
    slong i;

    if (vec->length > len)
    {
        for (i = len; i < vec->length; i++)
            GR_IGNORE(gr_poly_zero(vec->entries + i, ctx));
    }
    else if (vec->length < len)
    {
        gr_poly_vec_fit_length(vec, len, ctx);
        for (i = vec->length; i < len; i++)
            GR_IGNORE(gr_poly_zero(vec->entries + i, ctx));
    }

    vec->length = len;
}

int
gr_poly_vec_append(gr_poly_vec_t vec, const gr_poly_t f, gr_ctx_t ctx)
{
    gr_poly_vec_fit_length(vec, vec->length + 1, ctx);
    vec->length++;
    return gr_poly_set(vec->entries + vec->length - 1, f, ctx);
}
