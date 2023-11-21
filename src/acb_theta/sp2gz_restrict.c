/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
sp2gz_restrict(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    slong g1 = sp2gz_dim(res);
    fmpz_mat_t a, b, c, d;

    fmpz_mat_window_init(a, mat, 0, 0, g1, g1);
    fmpz_mat_window_init(b, mat, 0, g, g1, g + g1);
    fmpz_mat_window_init(c, mat, g, 0, g + g1, g1);
    fmpz_mat_window_init(d, mat, g, g, g + g1, g + g1);

    sp2gz_set_blocks(res, a, b, c, d);

    fmpz_mat_window_clear(a);
    fmpz_mat_window_clear(b);
    fmpz_mat_window_clear(c);
    fmpz_mat_window_clear(d);
}
