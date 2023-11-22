/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
sp2gz_is_correct(const fmpz_mat_t mat)
{
    slong r = fmpz_mat_nrows(mat);
    slong c = fmpz_mat_ncols(mat);
    slong g = r / 2;
    fmpz_mat_t J, test;
    int res;

    if (r != c || r % 2 != 0)
    {
        return 0;
    }

    fmpz_mat_init(J, 2 * g, 2 * g);
    fmpz_mat_init(test, 2 * g, 2 * g);

    sp2gz_j(J);
    fmpz_mat_transpose(test, mat);
    fmpz_mat_mul(test, test, J);
    fmpz_mat_mul(test, test, mat);

    res = fmpz_mat_equal(test, J);

    fmpz_mat_clear(J);
    fmpz_mat_clear(test);
    return res;
}
