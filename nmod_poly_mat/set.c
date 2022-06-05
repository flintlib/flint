/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly_mat-impl.h"


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
