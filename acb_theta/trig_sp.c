/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
fmpz_mat_trig_sp(fmpz_mat_t mat, const fmpz_mat_t S)
{
    slong g = fmpz_mat_nrows(mat) / 2;
    fmpz_mat_t zero, one;

    fmpz_mat_init(zero, g, g);
    fmpz_mat_init(one, g, g);

    fmpz_mat_one(one);
    fmpz_mat_set_abcd(mat, one, S, zero, one);

    fmpz_mat_clear(zero);
    fmpz_mat_clear(one);
}
