/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"

void
nmod_mat_concat_horizontal(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)
{
    slong i;
    slong r = mat1->r;
    slong c1 = mat1->c;
    slong c2 = mat2->c;

    for (i = 0; i < r; i++)
    {
    	flint_mpn_copyi(res->rows[i], mat1->rows[i], c1);
    	flint_mpn_copyi(res->rows[i] + c1, mat2->rows[i], c2);
    }
}

void
nmod_mat_concat_vertical(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)
{
    slong i;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong r2 = mat2->r;

    for (i = 0; i < r1; i++)
    	flint_mpn_copyi(res->rows[i], mat1->rows[i], c1);

    for (i = 0; i < r2; i++)
    	flint_mpn_copyi(res->rows[i + r1], mat2->rows[i], c1);
}
