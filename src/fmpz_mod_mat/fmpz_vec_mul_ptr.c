/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

void fmpz_mod_mat_fmpz_vec_mul_ptr(fmpz * const * c,
                   const fmpz * const * a, slong alen, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mat_fmpz_vec_mul_ptr(c, a, alen, B);
    for (i = 0; i < fmpz_mod_mat_ncols(B, ctx); i++)
        fmpz_mod_set_fmpz(c[i], c[i], ctx);
}
