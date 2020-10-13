/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"

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
