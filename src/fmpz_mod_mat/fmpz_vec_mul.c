/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"

void fmpz_mod_mat_fmpz_vec_mul(fmpz * c, const fmpz * a, slong alen,
                                                        const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_fmpz_vec_mul(c, a, alen, B);
    _fmpz_mod_vec_set_fmpz_vec(c, c, fmpz_mod_mat_ncols(B, ctx), ctx);
}
