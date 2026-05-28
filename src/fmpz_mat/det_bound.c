/*
    Copyright (C) 2011, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "mag.h"

static void
fmpz_mat_det_bound_inner(fmpz_t bound, const fmpz_mat_t A, int include_zero_norms)
{
    slong r = fmpz_mat_nrows(A);
    slong c = fmpz_mat_ncols(A);
    slong i, j;
    mag_t t, rbound, cbound;
    mag_struct *rnorms, *cnorms;

    if (r == 0 || c == 0)
    {
        fmpz_one(bound);
        return;
    }

    rnorms = _mag_vec_init(r + c);
    cnorms = rnorms + r;

    mag_init(t);
    mag_init(rbound);
    mag_init(cbound);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            mag_set_fmpz(t, fmpz_mat_entry(A, i, j));
            mag_fast_addmul(rnorms + i, t, t);
            mag_fast_addmul(cnorms + j, t, t);
        }
    }

    mag_one(rbound);
    mag_one(cbound);

    for (i = 0; i < r; i++)
    {
        if (include_zero_norms || !mag_is_zero(rnorms + i))
            mag_mul(rbound, rbound, rnorms + i);
    }

    for (i = 0; i < c; i++)
    {
        if (include_zero_norms || !mag_is_zero(cnorms + i))
            mag_mul(cbound, cbound, cnorms + i);
    }

    mag_min(t, rbound, cbound);
    mag_sqrt(t, t);
    mag_get_fmpz(bound, t);

    _mag_vec_clear(rnorms, r + c);

    mag_clear(t);
    mag_clear(rbound);
    mag_clear(cbound);
}


void
fmpz_mat_det_bound(fmpz_t bound, const fmpz_mat_t A)
{
    fmpz_mat_det_bound_inner(bound, A, 1);
}

void
fmpz_mat_det_bound_submatrix(fmpz_t bound, const fmpz_mat_t A)
{
    fmpz_mat_det_bound_inner(bound, A, 0);
}

