/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_randtest_spd(arb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits)
{
    slong g = arb_mat_nrows(mat);
    arb_mat_t tp;

    arb_mat_init(tp, g, g);
    arb_mat_randtest_cho(mat, state, prec, mag_bits);
    arb_mat_transpose(tp, mat);
    arb_mat_mul(mat, mat, tp, prec);

    arb_mat_clear(tp);
}
