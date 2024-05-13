/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010, 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod_mat.h"

void
nmod_mat_init(nmod_mat_t mat, slong rows, slong cols, ulong n)
{
    slong i;

    if (rows != 0)
        mat->rows = (ulong **) flint_malloc(rows * sizeof(ulong *));
    else
        mat->rows = NULL;

    if (rows != 0 && cols != 0)
    {
        mat->entries = (ulong *) flint_calloc(flint_mul_sizes(rows, cols), sizeof(ulong));

        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
    {
        mat->entries = NULL;
	if (rows != 0)
        {
            for (i = 0; i < rows; i++)
                mat->rows[i] = NULL;
        }
    }

    mat->r = rows;
    mat->c = cols;

    nmod_mat_set_mod(mat, n);
}

void
nmod_mat_init_set(nmod_mat_t mat, const nmod_mat_t src)
{
    slong rows = src->r;
    slong cols = src->c;
    slong i;

    if (rows != 0)
        mat->rows = flint_malloc(rows * sizeof(ulong *));
    else
        mat->rows = NULL;

    if ((rows) && (cols))
    {
        mat->entries = flint_malloc(flint_mul_sizes(rows, cols) * sizeof(ulong));

        for (i = 0; i < rows; i++)
        {
            mat->rows[i] = mat->entries + i * cols;
            flint_mpn_copyi(mat->rows[i], src->rows[i], cols);
        }
    }
    else
    {
        mat->entries = NULL;
	if (rows != 0)
        {
            for (i = 0; i < rows; i++)
                mat->rows[i] = NULL;
	}
    }

    mat->r = rows;
    mat->c = cols;

    mat->mod = src->mod;
}
