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
sp2gz_is_trig(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t alpha, gamma;
    int res;

    fmpz_mat_window_init(alpha, mat, 0, 0, g, g);
    fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);

    res = fmpz_mat_is_one(alpha) && fmpz_mat_is_zero(gamma);

    fmpz_mat_window_clear(alpha);
    fmpz_mat_window_clear(gamma);
    return res;
}
