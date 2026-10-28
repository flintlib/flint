/*
    Copyright (C) 2015 Dharak Kharod
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void fmpz_mat_content(fmpz_t mat_gcd, const fmpz_mat_t A)
{
    slong i;

    fmpz_zero(mat_gcd);

    for (i = 0; i < A->r && !fmpz_is_one(mat_gcd); i++)
        _fmpz_vec_content_chained(mat_gcd, fmpz_mat_row(A, i), A->c, mat_gcd);
}
