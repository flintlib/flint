/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void fmpz_mod_mat_scalar_mul_fmpz(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, fmpz_t c)
{
    fmpz_t d;
    fmpz_init(d);
    fmpz_mod(d, c, A->mod);

    fmpz_mat_scalar_mul_fmpz(B->mat, A->mat, d);
    _fmpz_mod_mat_reduce(B);
    fmpz_clear(d);
}
