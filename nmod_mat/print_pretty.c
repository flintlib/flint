/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

void
nmod_mat_print_pretty(const nmod_mat_t mat)
{
    slong i, j;
    int width;
    char fmt[FLINT_BITS + 5];

    flint_printf("<%wd x %wd integer matrix mod %wu>\n", mat->r, mat->c, mat->mod.n);

    if (!(mat->c) || !(mat->r))
        return;

    width = n_sizeinbase(mat->mod.n, 10);

    flint_sprintf(fmt, "%%%dwu", width);

    for (i = 0; i < mat->r; i++)
    {
        flint_printf("[");

        for (j = 0; j < mat->c; j++)
        {
            flint_printf(fmt, mat->rows[i][j]);
            if (j + 1 < mat->c)
                flint_printf(" ");
        }

        flint_printf("]\n");
    }
}

