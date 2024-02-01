/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void nmod_poly_mat_shift_left(nmod_poly_mat_t res, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_shift_left(nmod_poly_mat_entry(res, i, j), nmod_poly_mat_entry(pmat, i, j), k);
}

void nmod_poly_mat_shift_right(nmod_poly_mat_t res, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_shift_right(nmod_poly_mat_entry(res, i, j), nmod_poly_mat_entry(pmat, i, j), k);
}
