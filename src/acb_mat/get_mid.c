/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"

void acb_mat_get_mid(acb_mat_t B, const acb_mat_t A)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
    {
        for (j = 0; j < acb_mat_ncols(A); j++)
        {
            arb_get_mid_arb(acb_realref(acb_mat_entry(B, i, j)), acb_realref(acb_mat_entry(A, i, j)));
            arb_get_mid_arb(acb_imagref(acb_mat_entry(B, i, j)), acb_imagref(acb_mat_entry(A, i, j)));
        }
    }
}
