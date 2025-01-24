/*
    Copyright (C) 2017 Luca De Feo
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"

void fmpz_mod_mat_add(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    slong i;
    slong r = fmpz_mod_mat_nrows(A, ctx);
    slong c = fmpz_mod_mat_ncols(A, ctx);

    if (c != 0)
        for (i = 0; i < r; i++)
            _fmpz_mod_vec_add(fmpz_mod_mat_row(C, i), fmpz_mod_mat_row(A, i), fmpz_mod_mat_row(B, i), c, ctx);
}
