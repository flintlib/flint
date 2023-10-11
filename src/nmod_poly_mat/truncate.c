/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly_mat.h"

void nmod_poly_mat_set_trunc(nmod_poly_mat_t tmat, const nmod_poly_mat_t pmat, long len)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_set_trunc(tmat->rows[i] + j, pmat->rows[i] + j, len);
}

void nmod_poly_mat_truncate(nmod_poly_mat_t pmat, long len)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_truncate(pmat->rows[i] + j, len);
}

void nmod_poly_mat_shift_left(nmod_poly_mat_t smat, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < smat->r; i++)
        for (slong j = 0; j < smat->c; j++)
            nmod_poly_shift_left(nmod_poly_mat_entry(smat, i, j), nmod_poly_mat_entry(pmat, i, j), k);
}

void nmod_poly_mat_shift_right(nmod_poly_mat_t smat, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < smat->r; i++)
        for (slong j = 0; j < smat->c; j++)
            nmod_poly_shift_right(smat->rows[i] + j, pmat->rows[i] + j, k);
}

