/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_set_round_arb_mat(acb_mat_t dest, const arb_mat_t src, slong prec)
{
    slong i, j;

    if (acb_mat_ncols(dest) != 0)
    {
        for (i = 0; i < acb_mat_nrows(dest); i++)
            for (j = 0; j < acb_mat_ncols(dest); j++)
                acb_set_round_arb(acb_mat_entry(dest, i, j),
                    arb_mat_entry(src, i, j), prec);
    }
}
