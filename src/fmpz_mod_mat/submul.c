/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void fmpz_mod_mat_submul(fmpz_mod_mat_t D,
                         const fmpz_mod_mat_t C,
                         const fmpz_mod_mat_t A,
                         const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_mat_t tmp;
    fmpz_mod_mat_init(tmp, A->r, B->c, ctx);
    fmpz_mod_mat_mul(tmp, A, B, ctx);
    fmpz_mod_mat_sub(D, C, tmp, ctx);
    fmpz_mod_mat_clear(tmp, ctx);
}
