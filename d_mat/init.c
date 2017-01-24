/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "d_mat.h"

void
d_mat_init(d_mat_t mat, slong rows, slong cols)
{
    if (rows != 0 && cols != 0)       /* Allocate space for r*c small entries */
    {
        slong i;
        mat->entries = (double *) flint_calloc(rows*cols, sizeof(double));
        mat->rows = (double **) flint_malloc(rows*sizeof(double *));  /* Initialise rows */

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
}
