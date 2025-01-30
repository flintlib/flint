/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"

void
acb_mat_swap_rows(acb_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s)
    {
        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        _acb_vec_swap(acb_mat_entry(mat, r, 0), acb_mat_entry(mat, s, 0), mat->c);
    }
}