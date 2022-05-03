/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "fmpz.h"

void
nmod_mat_scalar_mul_fmpz(nmod_mat_t res, const nmod_mat_t M, const fmpz_t c)
{
    fmpz_t d;
    fmpz_init(d);
    fmpz_mod_ui(d, c, res->mod.n);
    nmod_mat_scalar_mul(res, M, fmpz_get_ui(d));
    fmpz_clear(d);
}
