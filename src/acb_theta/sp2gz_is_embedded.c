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
sp2gz_is_embedded(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t t;
    int r;

    fmpz_mat_init(t, 2 * g, 2 * g);

    sp2gz_restrict(res, mat);
    sp2gz_embed(t, res);
    r = fmpz_mat_equal(t, mat);

    fmpz_mat_clear(t);
    return r;
}
