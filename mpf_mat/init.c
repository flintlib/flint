/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "mpf_mat.h"

void
mpf_mat_init(mpf_mat_t mat, slong rows, slong cols, flint_bitcnt_t prec)
{

    if (rows != 0 && cols != 0)       /* Allocate space for r*c small entries */
    {
        slong i;
        mat->entries = (mpf *) flint_malloc(flint_mul_sizes(rows, cols) * sizeof(mpf));
        mat->rows = (mpf **) flint_malloc(rows * sizeof(mpf *));    /* Initialise rows */

        for (i = 0; i < rows * cols; i++)
            mpf_init2(mat->entries + i, prec);
        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
    {
       mat->entries = NULL;
       mat->rows = NULL;
    }

    mat->r = rows;
    mat->c = cols;
    mat->prec = prec;
}
