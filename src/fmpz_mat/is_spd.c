/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

int
fmpz_mat_is_spd(const fmpz_mat_t A)
{
    slong d = fmpz_mat_nrows(A);
    fmpz_mat_t tp, w;
    fmpz_t det;
    slong k;
    int res = 1;

    if (fmpz_mat_ncols(A) != d)
    {
        return 0;
    }

    fmpz_mat_init(tp, d, d);
    fmpz_init(det);

    fmpz_mat_transpose(tp, A);
    if (!fmpz_mat_equal(tp, A))
    {
        res = 0;
    }

    for (k = 1; (k <= d) && res; k++)
    {
        fmpz_mat_window_init(w, A, 0, 0, k, k);
        fmpz_mat_det(det, w);
        res = (fmpz_cmp_si(det, 0) > 0);
        fmpz_mat_window_clear(w);
    }

    fmpz_mat_clear(tp);
    fmpz_clear(det);
    return res;
}
