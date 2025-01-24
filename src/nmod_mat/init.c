/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010, 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "long_extras.h"
#include "mpn_extras.h"
#include "nmod_mat.h"

/* Set to 2 to check nonstandard strides */
#define STRIDE_DEBUG 1

void
nmod_mat_init(nmod_mat_t mat, slong rows, slong cols, ulong n)
{
    mat->r = rows;
    mat->c = cols;
    mat->entries = NULL;
    mat->stride = STRIDE_DEBUG * cols;

    if (rows != 0 && cols != 0)
    {
        slong num;
        int of;

        of = z_mul_checked(&num, rows, cols);

        if (of)
            flint_throw(FLINT_ERROR, "Overflow creating a %wd x %wd object\n", rows, cols);

        mat->entries = flint_calloc(STRIDE_DEBUG * num, sizeof(ulong));
    }

    nmod_mat_set_mod(mat, n);
}

void
nmod_mat_init_set(nmod_mat_t mat, const nmod_mat_t src)
{
    slong rows = src->r;
    slong cols = src->c;
    slong i;

    mat->r = rows;
    mat->c = cols;
    mat->entries = NULL;
    mat->stride = cols;
    mat->mod = src->mod;

    if (rows != 0 && cols != 0)
    {
        slong num;
        int of;

        of = z_mul_checked(&num, rows, cols);

        if (of)
            flint_throw(FLINT_ERROR, "Overflow creating a %wd x %wd object\n", rows, cols);

        mat->entries = flint_malloc(num * sizeof(ulong));

        for (i = 0; i < rows; i++)
            flint_mpn_copyi(nmod_mat_entry_ptr(mat, i, 0), nmod_mat_entry_ptr(src, i, 0), cols);
    }
}
