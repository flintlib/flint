/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "arb_mat.h"

void
arb_mat_entrywise_not_is_zero(fmpz_mat_t dest, const arb_mat_t src)
{
    slong i, j;
    fmpz_mat_zero(dest);
    for (i = 0; i < arb_mat_nrows(src); i++)
    {
        for (j = 0; j < arb_mat_ncols(src); j++)
        {
            if (!arb_is_zero(arb_mat_entry(src, i, j)))
            {
                fmpz_one(fmpz_mat_entry(dest, i, j));
            }
        }
    }
}
