/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_conjugate(acb_mat_t B, const acb_mat_t A)
{
    slong i, j;

    if ((acb_mat_nrows(B) != acb_mat_nrows(A)) ||
        (acb_mat_ncols(B) != acb_mat_ncols(A)))
    {
        flint_throw(FLINT_ERROR, "acb_mat_conjugate: incompatible dimensions.\n");
    }

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_conj(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j));
}

