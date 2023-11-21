/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_print(const ca_mat_t mat, ca_ctx_t ctx)
{
    slong r, c;
    slong i, j;

    r = ca_mat_nrows(mat);
    c = ca_mat_ncols(mat);

    flint_printf("ca_mat of size %wd x %wd:\n", r, c);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            flint_printf("    ");
            ca_print(ca_mat_entry(mat, i, j), ctx);
            flint_printf("\n");
        }

    }

    flint_printf("\n");
}

void
ca_mat_printn(const ca_mat_t mat, slong digits, ca_ctx_t ctx)
{
    slong r, c;
    slong i, j;

    r = ca_mat_nrows(mat);
    c = ca_mat_ncols(mat);

    flint_printf("[");

    for (i = 0; i < r; i++)
    {
        flint_printf("[");

        for (j = 0; j < c; j++)
        {
            ca_printn(ca_mat_entry(mat, i, j), digits, ctx);
            if (j < c - 1)
                flint_printf(", ");
        }

        if (i < r - 1)
            flint_printf("],\n");
        else
            flint_printf("]");
    }
    flint_printf("]\n");
}
