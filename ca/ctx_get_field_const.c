/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

slong
_ca_ctx_get_field_const(ca_ctx_t ctx, ulong func)
{
    slong i;

    for (i = 0; i < ctx->fields_len; i++)
    {
        /* todo: check length? (if symbols for constants are ever overloaded to take parameters...) */
        if (((ctx->fields + i)->type == CA_FIELD_TYPE_FUNC) && ((ctx->fields + i)->data.func.func == func))
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
        ca_field_init_const(ctx->fields + i, func);
    }

    return i;
}
