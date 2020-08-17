/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpf_mat.h"

void
mpf_mat_print(const mpf_mat_t mat)
{
    slong i, j;

    flint_printf("[");
    for (i = 0; i < mat->r; i++)
    {
        flint_printf("[");
        for (j = 0; j < mat->c; j++)
        {
            mpf_out_str(stdout, 10, 0, mpf_mat_entry(mat, i, j));
            if (j < mat->c - 1)
                flint_printf(" ");
        }
        flint_printf("]\n");
    }
    flint_printf("]\n");
}
