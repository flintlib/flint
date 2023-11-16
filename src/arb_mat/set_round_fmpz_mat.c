/*
    Copyright (C) 2012 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "arb_mat.h"

void
arb_mat_set_round_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src, slong prec)
{
    slong i, j;

    if (arb_mat_ncols(dest) != 0)
    {
        for (i = 0; i < arb_mat_nrows(dest); i++)
            for (j = 0; j < arb_mat_ncols(dest); j++)
                arb_set_round_fmpz(arb_mat_entry(dest, i, j),
                    fmpz_mat_entry(src, i, j), prec);
    }
}
