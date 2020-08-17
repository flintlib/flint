/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void fmpz_mod_mat_neg(fmpz_mod_mat_t B, const fmpz_mod_mat_t A)
{
    fmpz_mat_neg(B->mat, A->mat);
    _fmpz_mod_mat_reduce(B);
}
