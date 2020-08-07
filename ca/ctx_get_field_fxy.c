/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

slong
_ca_ctx_get_field_fxy(ca_ctx_t ctx, calcium_func_code func, const ca_t x, const ca_t y)
{
    slong i;

    for (i = 0; i < ctx->fields_len; i++)
    {
        if (CA_FIELD_LENGTH(ctx->fields + i) == 1 &&
            CA_EXT_HEAD(CA_FIELD_EXT_ELEM(ctx->fields + i, 0)) == func &&
            CA_EXT_FUNC_NARGS(CA_FIELD_EXT_ELEM(ctx->fields + i, 0)) == 2 &&
            ca_equal_repr(x, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(ctx->fields + i, 0)), ctx) &&
            ca_equal_repr(y, CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(ctx->fields + i, 0)) + 1, ctx))
        {
            break;
        }
    }

    if (i >= ctx->fields_len)
    {
        if (i >= ctx->fields_alloc)
        {
            ctx->fields = (ca_field_struct *) flint_realloc(ctx->fields, sizeof(ca_field_struct) * 2 * ctx->fields_alloc);
            ctx->fields_alloc = 2 * ctx->fields_alloc;
        }

        ctx->fields_len = i + 1;
        ca_field_init_fxy(ctx->fields + i, func, x, y, ctx);
    }

    return i;
}

