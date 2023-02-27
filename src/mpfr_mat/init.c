/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_mat.h"

void
mpfr_mat_init(mpfr_mat_t mat, slong rows, slong cols, mpfr_prec_t prec)
{

    if (rows != 0 && cols != 0)       /* Allocate space for r*c small entries */
    {
        slong i;
        mat->entries =
            (__mpfr_struct *) flint_malloc(flint_mul_sizes(rows, cols) * sizeof(__mpfr_struct));
        mat->rows = (__mpfr_struct **) flint_malloc(rows * sizeof(__mpfr_struct *));  /* Initialise rows */

        for (i = 0; i < rows * cols; i++)
            mpfr_init2(mat->entries + i, prec);
        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
        mat->entries = NULL;

    mat->r = rows;
    mat->c = cols;
    mat->prec = prec;
}
