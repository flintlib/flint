/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_field_print(const ca_field_t K, const ca_ctx_t ctx)
{
    slong i, len;

    flint_printf("QQ");

    if (K->type == CA_FIELD_TYPE_QQ)
        return;

    if (K->type == CA_FIELD_TYPE_NF)
    {
        flint_printf("(x) where {x = Algebraic number ");
        qqbar_printnd(&K->data.nf.x, 10);
        flint_printf("}");
        return;
    }

    if (K->type == CA_FIELD_TYPE_FUNC)
    {
        flint_printf("(x) where {x = ");
        flint_printf("%s", calcium_func_name(K->data.func.func));
        if (K->data.func.args_len != 0)
        {
            slong i;
            flint_printf("(");
            for (i = 0; i < K->data.func.args_len; i++)
            {
                ca_print(K->data.func.args + i, ctx);
                if (i < K->data.func.args_len - 1)
                    flint_printf(", ");
            }
            flint_printf(")");
        }
        flint_printf("}");
    }

    if (K->type == CA_FIELD_TYPE_MULTI)
    {
        len = K->data.multi.len;

        flint_printf("(");
        for (i = 0; i < len; i++)
        {
            flint_printf("x%wd", i + 1);
            if (i < len - 1)
                flint_printf(", ");
        }
        flint_printf(") where {");
        for (i = 0; i < len; i++)
        {
            flint_printf("x%wd = [", i + 1);

            ca_field_print(ctx->fields + K->data.multi.ext[i], ctx);

            flint_printf("]");

            if (i < len - 1)
                flint_printf(", ");
        }
        flint_printf("}");

        if (K->data.multi.ideal_len > 0)
        {
            /* todo: could use cached mctx in ctx... */
            fmpz_mpoly_ctx_t mctx;
            fmpz_mpoly_ctx_init(mctx, K->data.multi.len, CA_MPOLY_ORD);

            flint_printf(" with ideal {");
            for (i = 0; i < K->data.multi.ideal_len; i++)
            {
                fmpz_mpoly_print_pretty(K->data.multi.ideal + i, NULL, mctx);
                if (i < K->data.multi.ideal_len - 1)
                    flint_printf(", ");
            }
            flint_printf("}");
        }
    }
}

