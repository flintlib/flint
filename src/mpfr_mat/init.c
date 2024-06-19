/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "long_extras.h"
#include "mpfr_mat.h"

void
mpfr_mat_init(mpfr_mat_t mat, slong rows, slong cols, mpfr_prec_t prec)
{
    mat->r = rows;
    mat->c = cols;
    mat->prec = prec;

    if (rows != 0 && cols != 0)
    {
        slong i;
        slong num;
        int of;

        of = z_mul_checked(&num, rows, cols);

        if (of)
            flint_throw(FLINT_ERROR, "Overflow creating a %wd x %wd object\n", rows, cols);

        mat->entries = flint_malloc(num * sizeof(__mpfr_struct));
        mat->rows = flint_malloc(rows * sizeof(__mpfr_struct *));

        for (i = 0; i < rows * cols; i++)
            mpfr_init2(mat->entries + i, prec);
        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
        mat->entries = NULL;
}
