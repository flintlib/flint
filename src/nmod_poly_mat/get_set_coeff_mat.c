/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void nmod_poly_mat_get_coeff_mat(nmod_mat_t res, const nmod_poly_mat_t mat, slong deg)
{
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_mat_set_entry(res, i, j,
                    nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j), deg));
}

void nmod_poly_mat_set_coeff_mat(nmod_poly_mat_t pmat,
                                 const nmod_mat_t coeff,
                                 slong deg)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(pmat, i, j),
                    deg, nmod_mat_entry(coeff, i, j));
}
