/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

static truth_t
ca_field_is_algebraic(const ca_field_t K, ca_ctx_t ctx)
{
    if (K->type == CA_FIELD_TYPE_NF)
        return T_TRUE;

    return T_UNKNOWN;
}

truth_t
ca_check_is_algebraic(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        if (ca_is_unknown(x, ctx))
            return T_UNKNOWN;

        return T_FALSE;
    }

    if (x->field == CA_FIELD_ID_QQ ||
        ctx->fields[x->field].type == CA_FIELD_TYPE_NF)
        return T_TRUE;

    if (ctx->fields[x->field].type == CA_FIELD_TYPE_MULTI)
    {
        slong len, i;

        len = ctx->fields[x->field].data.multi.len;

        /* todo: handle simple transcendental numbers, e.g. Q(i,pi) */
        /* need to verify that some the generator is used in the poly */
        /* for Q(a,b,pi) we don't know, because a, b could cancel out pi */
        for (i = 0; len; i++)
        {
            if (ca_field_is_algebraic(ctx->fields + ctx->fields[x->field].data.multi.ext[i], ctx) != T_TRUE)
                return T_UNKNOWN;
        }

        return T_TRUE;
    }

    return T_UNKNOWN;
}

