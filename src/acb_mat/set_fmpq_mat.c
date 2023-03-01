/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_set_fmpq_mat(acb_mat_t dest, const fmpq_mat_t src, slong prec)
{
    slong i, j;

    if (acb_mat_ncols(dest) != 0)
    {
        for (i = 0; i < acb_mat_nrows(dest); i++)
            for (j = 0; j < acb_mat_ncols(dest); j++)
                acb_set_fmpq(acb_mat_entry(dest, i, j),
                    fmpq_mat_entry(src, i, j), prec);
    }
}
