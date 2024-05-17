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

void acb_mat_add_error_mag(acb_mat_t mat, const mag_t err)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(mat); i++)
        for (j = 0; j < acb_mat_ncols(mat); j++)
            acb_add_error_mag(acb_mat_entry(mat, i, j), err);
}
