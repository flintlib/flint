/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"
#include "arb_mat.h"

void
arb_mat_set_fmpq_mat(arb_mat_t dest, const fmpq_mat_t src, slong prec)
{
    slong i, j;

    if (arb_mat_ncols(dest) != 0)
    {
        for (i = 0; i < arb_mat_nrows(dest); i++)
            for (j = 0; j < arb_mat_ncols(dest); j++)
                arb_set_fmpq(arb_mat_entry(dest, i, j),
                    fmpq_mat_entry(src, i, j), prec);
    }
}
