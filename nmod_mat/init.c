/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"


void
nmod_mat_init(nmod_mat_t mat, slong rows, slong cols, mp_limb_t n)
{
    if (rows > 0 && cols > 0)
    {
        slong i;
        mat->entries = (mp_limb_t *) flint_calloc(flint_mul_sizes(rows, cols), sizeof(mp_limb_t));
        mat->rows = (mp_limb_t **) flint_malloc(rows * sizeof(mp_limb_t *));

        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;

        mat->r = rows;
        mat->c = cols;
    }
    else
        memset(mat, 0, sizeof(*mat));

    _nmod_mat_set_mod(mat, n);
}
