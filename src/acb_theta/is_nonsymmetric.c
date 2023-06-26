/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
arb_mat_is_nonsymmetric(const arb_mat_t mat)
{
    arb_mat_t tp;
    slong nrows = arb_mat_nrows(mat);
    int res;

    if (nrows != arb_mat_ncols(mat))
        return 1;

    arb_mat_init(tp, nrows, nrows);
    arb_mat_transpose(tp, mat);
    res = !arb_mat_overlaps(tp, mat);
    arb_mat_clear(tp);

    return res;
}
