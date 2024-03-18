/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_mat.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_set(nmod_poly_mat_t B, const nmod_poly_mat_t A)
{
    if (A != B)
    {
        slong i, j;

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                nmod_poly_set(nmod_poly_mat_entry(B, i, j),
                              nmod_poly_mat_entry(A, i, j));
    }
}

void nmod_poly_mat_set_nmod_mat(nmod_poly_mat_t pmat, const nmod_mat_t cmat)
{
    for (slong i = 0; i < cmat->r; i++)
    {
        for (slong j = 0; j < cmat->c; j++)
        {
            if (nmod_mat_entry(cmat, i, j) == 0)
                nmod_poly_zero(nmod_poly_mat_entry(pmat, i, j));
            else
            {
                nmod_poly_realloc(nmod_poly_mat_entry(pmat, i, j), 1);
                nmod_poly_mat_entry(pmat, i, j)->coeffs[0]
                    = nmod_mat_entry(cmat, i, j);
            }
        }
    }
}
