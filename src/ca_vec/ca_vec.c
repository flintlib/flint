/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"

void
_ca_vec_add(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_add(res + i, vec1 + i, vec2 + i, ctx);
}

truth_t
_ca_vec_check_is_zero(ca_srcptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    int have_unknown;
    truth_t is_zero;

    have_unknown = 0;
    for (i = 0; i < len; i++)
    {
        is_zero = ca_check_is_zero(vec + i, ctx);

        if (is_zero == T_FALSE)
            return T_FALSE;

        if (is_zero == T_UNKNOWN)
            have_unknown = 1;
    }

    if (have_unknown)
        return T_UNKNOWN;
    else
        return T_TRUE;
}

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

void
_ca_vec_neg(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_neg(res + i, src + i, ctx);
}

void
ca_vec_neg(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx)
{
    if (res != src)
    {
        ca_vec_set_length(res, ca_vec_length(src, ctx), ctx);
        _ca_vec_neg(ca_vec_entry(res, 0), ca_vec_entry(src, 0), ca_vec_length(res, ctx), ctx);
    }
}

void
ca_vec_print(const ca_vec_t vec, ca_ctx_t ctx)
{
    slong i, len;

    len = vec->length;

    flint_printf("ca_vec of length %wd:\n", len);

    for (i = 0; i < len; i++)
    {
        flint_printf("    ");
        ca_print(vec->entries + i, ctx);
        flint_printf("\n");
    }

    flint_printf("\n");
}

void
ca_vec_printn(const ca_vec_t vec, slong digits, ca_ctx_t ctx)
{
    slong len, i;

    len = vec->length;

    flint_printf("[");

    for (i = 0; i < len; i++)
    {
        ca_printn(vec->entries + i, digits, ctx);
        if (i < len - 1)
            flint_printf(", ");
    }
    flint_printf("]\n");
}

void
_ca_vec_scalar_addmul_ca(ca_ptr res, ca_srcptr vec, slong len, const ca_t c, ca_ctx_t ctx)
{
    slong i;
    ca_t t;

    if (len > 0)
    {
        ca_init(t, ctx);
        for (i = 0; i < len; i++)
        {
            ca_mul(t, vec + i, c, ctx);
            ca_add(res + i, res + i, t, ctx);
        }
        ca_clear(t, ctx);
    }
}

void
_ca_vec_scalar_div_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)
{
    slong i;

    if (len > 0)
    {
        ca_t t;
        ca_init(t, ctx);
        ca_inv(t, c, ctx);

        for (i = 0; i < len; i++)
            ca_mul(res + i, src + i, t, ctx);
        ca_clear(t, ctx);
    }
}

void
_ca_vec_scalar_mul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_mul(res + i, src + i, c, ctx);
}

void
_ca_vec_scalar_submul_ca(ca_ptr res, ca_srcptr vec, slong len, const ca_t c, ca_ctx_t ctx)
{
    slong i;
    ca_t t;

    if (len > 0)
    {
        ca_init(t, ctx);
        for (i = 0; i < len; i++)
        {
            ca_mul(t, vec + i, c, ctx);
            ca_sub(res + i, res + i, t, ctx);
        }
        ca_clear(t, ctx);
    }
}

void
_ca_vec_set(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_set(res + i, src + i, ctx);
}

void
ca_vec_set(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx)
{
    if (res != src)
    {
        ca_vec_set_length(res, ca_vec_length(src, ctx), ctx);
        _ca_vec_set(ca_vec_entry(res, 0), ca_vec_entry(src, 0), ca_vec_length(res, ctx), ctx);
    }
}

void
_ca_vec_fit_length(ca_vec_t vec, slong len, ca_ctx_t ctx)
{
    if (len > vec->alloc)
    {
        slong i;

        if (len < 2 * vec->alloc)
            len = 2 * vec->alloc;

        vec->entries = flint_realloc(vec->entries, len * sizeof(ca_struct));

        for (i = vec->alloc; i < len; i++)
            ca_init(ca_vec_entry(vec, i), ctx);

        vec->alloc = len;
    }
}

void
ca_vec_set_length(ca_vec_t vec, slong len, ca_ctx_t ctx)
{
    slong i;

    if (vec->length > len)
    {
        for (i = len; i < vec->length; i++)
            ca_zero(ca_vec_entry(vec, i), ctx);
    }
    else if (vec->length < len)
    {
        _ca_vec_fit_length(vec, len, ctx);

        for (i = vec->length; i < len; i++)
            ca_zero(ca_vec_entry(vec, i), ctx);
    }

    vec->length = len;
}

void
_ca_vec_sub(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_sub(res + i, vec1 + i, vec2 + i, ctx);
}

void
_ca_vec_swap(ca_ptr vec1, ca_ptr vec2, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_swap(vec1 + i, vec2 + i, ctx);
}

void
_ca_vec_zero(ca_ptr res, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        ca_zero(res + i, ctx);
}

void
ca_vec_zero(ca_vec_t res, slong len, ca_ctx_t ctx)
{
    ca_vec_set_length(res, len, ctx);
    _ca_vec_zero(ca_vec_entry(res, 0), ca_vec_length(res, ctx), ctx);
}

