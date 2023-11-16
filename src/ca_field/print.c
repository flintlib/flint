/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

void
ca_field_print(const ca_field_t K, ca_ctx_t ctx)
{
    slong i, len, ideal_len;

    flint_printf("QQ");

    len = CA_FIELD_LENGTH(K);

    if (len == 0)
        return;

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
        flint_printf("x%wd = ", i + 1);

        ca_ext_print(CA_FIELD_EXT_ELEM(K, i), ctx);

        flint_printf("");

        if (i < len - 1)
            flint_printf(", ");
    }
    flint_printf("}");

    ideal_len = CA_FIELD_IDEAL_LENGTH(K);

    if (ideal_len > 0)
    {
        flint_printf(" with ideal {");
        for (i = 0; i < ideal_len; i++)
        {
            fmpz_mpoly_print_pretty(CA_FIELD_IDEAL_ELEM(K, i), NULL, CA_FIELD_MCTX(K, ctx));
            if (i < ideal_len - 1)
                flint_printf(", ");
        }
        flint_printf("}");
    }
}

