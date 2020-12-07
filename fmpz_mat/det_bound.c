/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

static void
fmpz_mat_det_bound_inner(fmpz_t bound, const fmpz_mat_t A, int zero_cols)
{
    fmpz_t p, s, t;
    slong i, j;

    fmpz_init(p);
    fmpz_init(s);
    fmpz_init(t);
    fmpz_one(p);

    for (i = 0; i < A->r; i++)
    {
        fmpz_zero(s);

        for (j = 0; j < A->c; j++)
            fmpz_addmul(s, A->rows[i] + j, A->rows[i] + j);

        fmpz_sqrtrem(s, t, s);
        if (!fmpz_is_zero(t))
            fmpz_add_ui(s, s, UWORD(1));

        if (zero_cols || !fmpz_is_zero(s))
            fmpz_mul(p, p, s);
    }

    fmpz_set(bound, p);
    fmpz_clear(p);
    fmpz_clear(s);
    fmpz_clear(t);
}

void
fmpz_mat_det_bound(fmpz_t bound, const fmpz_mat_t A)
{
    fmpz_mat_det_bound_inner(bound, A, 1);
}

void
fmpz_mat_det_bound_nonzero(fmpz_t bound, const fmpz_mat_t A)
{
    fmpz_mat_det_bound_inner(bound, A, 0);
}
