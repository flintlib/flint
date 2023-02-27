/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_pascal(arb_mat_t mat, int triangular, slong prec)
{
    slong R, C, i, j;

    R = arb_mat_nrows(mat);
    C = arb_mat_ncols(mat);

    if (R == 0 || C == 0)
        return;

    if (triangular == 1)
    {
        for (i = 1; i < R; i++)
            for (j = 0; j < i && j < C; j++)
                arb_zero(arb_mat_entry(mat, i, j));

        for (j = 0; j < C; j++)
            arb_one(arb_mat_entry(mat, 0, j));

        for (i = 1; i < R && i < C; i++)
            arb_one(arb_mat_entry(mat, i, i));

        for (i = 1; i < R; i++)
            for (j = i + 1; j < C; j++)
                arb_add(arb_mat_entry(mat, i, j),
                    arb_mat_entry(mat, i, j - 1),
                    arb_mat_entry(mat, i - 1, j - 1), prec);
    }
    else if (triangular == -1)
    {
        for (i = 0; i < R; i++)
            for (j = i + 1; j < C; j++)
                arb_zero(arb_mat_entry(mat, i, j));

        for (i = 0; i < R; i++)
            arb_one(arb_mat_entry(mat, i, 0));

        for (i = 1; i < R && i < C; i++)
            arb_one(arb_mat_entry(mat, i, i));

        for (i = 2; i < R; i++)
            for (j = 1; j < i && j < C; j++)
                arb_add(arb_mat_entry(mat, i, j),
                    arb_mat_entry(mat, i - 1, j - 1),
                    arb_mat_entry(mat, i - 1, j), prec);
    }
    else
    {
        for (j = 0; j < C; j++)
            arb_one(arb_mat_entry(mat, 0, j));

        for (i = 1; i < R; i++)
            arb_one(arb_mat_entry(mat, i, 0));

        for (i = 1; i < R; i++)
            for (j = 1; j < C; j++)
                arb_add(arb_mat_entry(mat, i, j),
                    arb_mat_entry(mat, i, j - 1),
                    arb_mat_entry(mat, i - 1, j), prec);
    }
}

