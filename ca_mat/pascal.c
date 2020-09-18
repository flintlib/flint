/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_pascal(ca_mat_t mat, int triangular, ca_ctx_t ctx)
{
    slong R, C, i, j;

    R = ca_mat_nrows(mat);
    C = ca_mat_ncols(mat);

    if (R == 0 || C == 0)
        return;

    if (triangular == 1)
    {
        for (i = 1; i < R; i++)
            for (j = 0; j < i && j < C; j++)
                ca_zero(ca_mat_entry(mat, i, j), ctx);

        for (j = 0; j < C; j++)
            ca_one(ca_mat_entry(mat, 0, j), ctx);

        for (i = 1; i < R && i < C; i++)
            ca_one(ca_mat_entry(mat, i, i), ctx);

        for (i = 1; i < R; i++)
            for (j = i + 1; j < C; j++)
                ca_add(ca_mat_entry(mat, i, j),
                    ca_mat_entry(mat, i, j - 1),
                    ca_mat_entry(mat, i - 1, j - 1), ctx);
    }
    else if (triangular == -1)
    {
        for (i = 0; i < R; i++)
            for (j = i + 1; j < C; j++)
                ca_zero(ca_mat_entry(mat, i, j), ctx);

        for (i = 0; i < R; i++)
            ca_one(ca_mat_entry(mat, i, 0), ctx);

        for (i = 1; i < R && i < C; i++)
            ca_one(ca_mat_entry(mat, i, i), ctx);

        for (i = 2; i < R; i++)
            for (j = 1; j < i && j < C; j++)
                ca_add(ca_mat_entry(mat, i, j),
                    ca_mat_entry(mat, i - 1, j - 1),
                    ca_mat_entry(mat, i - 1, j), ctx);
    }
    else
    {
        for (j = 0; j < C; j++)
            ca_one(ca_mat_entry(mat, 0, j), ctx);

        for (i = 1; i < R; i++)
            ca_one(ca_mat_entry(mat, i, 0), ctx);

        for (i = 1; i < R; i++)
            for (j = 1; j < C; j++)
                ca_add(ca_mat_entry(mat, i, j),
                    ca_mat_entry(mat, i, j - 1),
                    ca_mat_entry(mat, i - 1, j), ctx);
    }
}
