/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2015 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_fprintd(FILE * file, const acb_mat_t mat, slong digits)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(mat); i++)
    {
        flint_fprintf(file, "[");

        for (j = 0; j < acb_mat_ncols(mat); j++)
        {
            acb_fprintd(file, acb_mat_entry(mat, i, j), digits);

            if (j < acb_mat_ncols(mat) - 1)
                flint_fprintf(file, ", ");
        }

        flint_fprintf(file, "]\n");
    }
}
