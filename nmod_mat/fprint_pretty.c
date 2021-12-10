/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include "nmod_mat.h"

int nmod_mat_fprint_pretty(FILE* file, const nmod_mat_t mat)
{
    slong i, j;
    int z, width;
    char fmt[FLINT_BITS + 5];

    z = flint_fprintf(file, "<%wd x %wd integer matrix mod %wu>\n",
                                                   mat->r, mat->c, mat->mod.n);
    if (z <= 0)
        return z;

    if (!(mat->c) || !(mat->r))
        return z;

    width = n_sizeinbase(mat->mod.n, 10);

    z = flint_sprintf(fmt, "%%%dwu", width);
    if (z <= 0)
        return z;

    for (i = 0; i < mat->r; i++)
    {
        z = flint_printf("[");
        if (z <= 0)
            return z;

        for (j = 0; j < mat->c; j++)
        {
            z = flint_printf(fmt, mat->rows[i][j]);
            if (z <= 0)
                return z;

            if (j + 1 < mat->c)
            {
                z = flint_printf(" ");
                if (z <= 0)
                    return z;
            }
        }

        flint_printf("]\n");
        if (z <= 0)
            return z;
    }

    return z;
}

