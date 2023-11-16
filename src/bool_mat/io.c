/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "bool_mat.h"

/* printing *******************************************************************/

void
bool_mat_fprint(FILE * file, const bool_mat_t mat)
{
    slong i, j;

    for (i = 0; i < bool_mat_nrows(mat); i++)
    {
        flint_fprintf(file, "[");

        for (j = 0; j < bool_mat_ncols(mat); j++)
        {
            flint_fprintf(file, "%d", bool_mat_get_entry(mat, i, j));

            if (j < bool_mat_ncols(mat) - 1)
                flint_fprintf(file, ", ");
        }

        flint_fprintf(file, "]\n");
    }
}

void bool_mat_print(const bool_mat_t mat) { bool_mat_fprint(stdout, mat); }
