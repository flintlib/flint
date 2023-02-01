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
sp2gz_is_correct(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t a, b, c, d;
    fmpz_mat_t prod1, prod2;
    int res;

    fmpz_mat_init(a, g, g);
    fmpz_mat_init(b, g, g);
    fmpz_mat_init(c, g, g);
    fmpz_mat_init(d, g, g);
    fmpz_mat_init(prod1, g, g);
    fmpz_mat_init(prod2, g, g);

    sp2gz_get_a(a, mat);
    sp2gz_get_b(b, mat);
    sp2gz_get_c(c, mat);
    sp2gz_get_d(d, mat);

    fmpz_mat_transpose(prod1, a);
    fmpz_mat_mul(prod1, prod1, c);
    fmpz_mat_transpose(prod2, c);
    fmpz_mat_mul(prod2, prod2, a);
    fmpz_mat_sub(prod1, prod1, prod2);
    res = fmpz_mat_is_zero(prod1);

    fmpz_mat_transpose(prod1, b);
    fmpz_mat_mul(prod1, prod1, d);
    fmpz_mat_transpose(prod2, d);
    fmpz_mat_mul(prod2, prod2, b);
    fmpz_mat_sub(prod1, prod1, prod2);
    res = res && fmpz_mat_is_zero(prod1);

    fmpz_mat_transpose(prod1, a);
    fmpz_mat_mul(prod1, prod1, d);
    fmpz_mat_transpose(prod2, c);
    fmpz_mat_mul(prod2, prod2, b);
    fmpz_mat_sub(prod1, prod1, prod2);
    res = res && fmpz_mat_is_one(prod1);

    fmpz_mat_clear(a);
    fmpz_mat_clear(b);
    fmpz_mat_clear(c);
    fmpz_mat_clear(d);
    fmpz_mat_clear(prod1);
    fmpz_mat_clear(prod2);

    return res;
}
