/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

void
fmpz_sparse_mat_det_bound(fmpz_t bound, const fmpz_sparse_mat_t A)
{
    fmpz_t s, t;
    slong i, j;

    fmpz_init(s);
    fmpz_init(t);

    /* bound = II_i ceil(||A[i]||_2) */
    fmpz_one(bound);
    for (i = 0; i < A->r; i++)
    {
        fmpz_zero(s);
        for (j = 0; j < A->rows[i].nnz; j++)
            fmpz_addmul(s, A->rows[i].entries[j].val, A->rows[i].entries[j].val);
        fmpz_sqrtrem(s, t, s);

        if (!fmpz_is_zero(t))
            fmpz_add_ui(s, s, UWORD(1));

        fmpz_mul(bound, bound, s);
    }

    fmpz_clear(s);
    fmpz_clear(t);
}
