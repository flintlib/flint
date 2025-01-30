/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mat.h"

void fmpz_mat_swap_rows(fmpz_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s && !fmpz_mat_is_empty(mat))
    {
        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        _fmpz_vec_swap(fmpz_mat_entry(mat, r, 0), fmpz_mat_entry(mat, s, 0), mat->c);
    }
}