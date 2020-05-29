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
ca_field_print(const ca_field_t K)
{
    slong i, len;
    len = K->len;

    flint_printf("QQ");
    if (len > 0)
    {
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
            if (K->type == CA_FIELD_TYPE_NF)
                ca_extension_print(K->nf_ext);
            else
                ca_extension_print(K->ext[i]);
            if (i < len - 1)
                flint_printf(", ");
        }
        flint_printf("}");

        if (K->ideal_len > 0)
        {
            flint_printf(" with ideal {");
            for (i = 0; i < K->ideal_len; i++)
            {
                fmpz_mpoly_print_pretty(K->ideal + i, NULL, CA_FIELD_MCTX(K));
                if (i < K->ideal_len - 1)
                    flint_printf(", ");
            }
            flint_printf("}");
        }
    }
}

